"""Local web UI: paste an ELP sequence, simulate it, read off the Rg.

Started with `sim ui`. Listens on localhost only. Jobs run one at a time in a
separate process (tools.rg_worker), each in its own folder under runs/.
"""

from __future__ import annotations

import csv
import json
import re
import secrets
import subprocess
import sys
import threading
from collections import deque
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

from tools import rg as rg_tools

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = PROJECT_ROOT / "runs"
PAGE = Path(__file__).with_name("ui.html")
JOB_ID = re.compile(r"\d{8}-\d{6}(-[a-z0-9-]+)?")
LIMITS = {"temp_c": (0.0, 100.0), "ionic": (0.01, 2.0), "ph": (1.0, 14.0)}
ERROR_TAIL_LINES = 25


class BadRequest(ValueError):
    pass


def validate_job(body: dict) -> dict:
    try:
        sequence = rg_tools.parse_sequence(str(body.get("sequence", "")))
    except rg_tools.SequenceError as exc:
        raise BadRequest(str(exc)) from exc

    preset = body.get("preset")
    if preset not in rg_tools.PRESETS:
        raise BadRequest("Choose a run length.")

    numbers = {}
    for key, (low, high) in LIMITS.items():
        try:
            value = float(body[key])
        except (KeyError, TypeError, ValueError):
            raise BadRequest(f"{key} must be a number.") from None
        if not low <= value <= high:
            raise BadRequest(f"{key} must be between {low:g} and {high:g}.")
        numbers[key] = value

    return {
        "name": str(body.get("name") or "").strip()[:80],
        "input": str(body["sequence"]).strip(),
        "sequence": sequence,
        "preset": preset,
        "plan": rg_tools.run_plan(len(sequence), preset),
        "seed": secrets.randbelow(2**31 - 1) + 1,
        **numbers,
    }


def read_progress(job_dir: Path, total_steps: int) -> dict:
    """Progress from the tab-separated StateDataReporter log: step, speed, elapsed."""
    try:
        lines = (job_dir / "elp.log").read_text().strip().splitlines()
        step, _, elapsed = lines[-1].split("\t")[:3]
        step, elapsed = int(step), float(elapsed)
    except (OSError, IndexError, ValueError):
        return {"fraction": 0.0, "eta_s": None}
    speed = step / elapsed if elapsed > 0 else 0
    eta = (total_steps - step) / speed if speed > 0 else None
    return {"fraction": min(step / total_steps, 1.0), "eta_s": eta}


class JobManager:
    def __init__(self, runs_dir: Path = RUNS_DIR) -> None:
        self.runs_dir = runs_dir
        self.runs_dir.mkdir(exist_ok=True)
        self.changed = threading.Condition()
        self.queued: deque[str] = deque()
        self.active: str | None = None
        self.proc: subprocess.Popen | None = None
        threading.Thread(target=self._loop, daemon=True).start()

    def submit(self, job: dict) -> str:
        slug = re.sub(r"[^a-z0-9]+", "-", job["name"].lower()).strip("-")[:40]
        stem = datetime.now().strftime("%Y%m%d-%H%M%S")
        job_id = f"{stem}-{slug}" if slug else stem
        job_dir = self.runs_dir / job_id
        suffix = 1
        while True:
            try:
                job_dir.mkdir()
                break
            except FileExistsError:
                suffix += 1
                job_id = f"{stem}-{slug}-{suffix}" if slug else f"{stem}-{suffix}"
                job_dir = self.runs_dir / job_id
        (job_dir / "job.json").write_text(json.dumps(job, indent=2))
        with self.changed:
            self.queued.append(job_id)
            self.changed.notify_all()
        return job_id

    def cancel(self, job_id: str) -> None:
        with self.changed:
            if job_id in self.queued:
                self.queued.remove(job_id)
                (self.runs_dir / job_id / "cancelled").touch()
            elif job_id == self.active and self.proc:
                (self.runs_dir / job_id / "cancelled").touch()
                self.proc.terminate()

    def shutdown(self) -> None:
        with self.changed:
            if self.active and self.proc:
                (self.runs_dir / self.active / "cancelled").touch()
                self.proc.terminate()

    def _loop(self) -> None:
        while True:
            with self.changed:
                while not self.queued:
                    self.changed.wait()
                self.active = job_id = self.queued.popleft()
                job_dir = self.runs_dir / job_id
                log = open(job_dir / "run.log", "wb")
                self.proc = subprocess.Popen(
                    [sys.executable, "-u", "-m", "tools.rg_worker", str(job_dir)],
                    cwd=PROJECT_ROOT,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
                proc = self.proc
            proc.wait()
            log.close()
            with self.changed:
                self.active, self.proc = None, None
                self._record_failure(job_dir)

    @staticmethod
    def _record_failure(job_dir: Path) -> None:
        if (job_dir / "result.json").exists() or (job_dir / "cancelled").exists():
            return
        tail = (job_dir / "run.log").read_text(errors="replace").splitlines()[-ERROR_TAIL_LINES:]
        (job_dir / "error.txt").write_text("\n".join(tail))

    def summary(self, job_id: str) -> dict:
        job_dir = self.runs_dir / job_id
        job = json.loads((job_dir / "job.json").read_text())
        out = {
            "id": job_id,
            "name": job["name"],
            "n_residues": len(job["sequence"]),
            "temp_c": job["temp_c"],
            "ionic": job["ionic"],
            "ph": job["ph"],
            "preset": job["preset"],
            "simulated_ns": job["plan"]["simulated_ns"],
        }
        with self.changed:
            active, queued = self.active, list(self.queued)
        if (job_dir / "result.json").exists():
            out["status"] = "done"
            out["result"] = json.loads((job_dir / "result.json").read_text())
        elif job_id == active:
            out["status"] = "running"
            out["progress"] = read_progress(job_dir, job["plan"]["steps"])
        elif job_id in queued:
            out["status"] = "queued"
        elif (job_dir / "cancelled").exists():
            out["status"] = "cancelled"
        elif (job_dir / "error.txt").exists():
            out["status"] = "failed"
            out["error"] = (job_dir / "error.txt").read_text()
        else:
            out["status"] = "interrupted"
        return out

    def list_jobs(self) -> list[dict]:
        ids = [
            d.name
            for d in self.runs_dir.iterdir()
            if d.is_dir() and JOB_ID.fullmatch(d.name) and (d / "job.json").exists()
        ]
        return [self.summary(i) for i in sorted(ids, reverse=True)]

    def detail(self, job_id: str) -> dict:
        out = self.summary(job_id)
        job_dir = self.runs_dir / job_id
        job = json.loads((job_dir / "job.json").read_text())
        out["sequence"] = job["sequence"]
        out["input"] = job["input"]
        if out["status"] == "done":
            with open(job_dir / "rg.csv", newline="") as f:
                rows = list(csv.DictReader(f))
            rg = np.array([float(r["rg_nm"]) for r in rows])
            skip = out["result"]["frames_discarded"]
            out["series"] = {
                "time_ns": [float(r["time_ns"]) for r in rows],
                "rg_nm": rg.tolist(),
                "cumulative_nm": [None] * skip + rg_tools.cumulative_mean(rg, skip)[skip:].tolist(),
            }
        return out


def make_handler(manager: JobManager, allowed_hosts: set[str]):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args) -> None:
            pass

        def _send(self, status: int, body: bytes, content_type: str, extra=None) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            for key, value in (extra or {}).items():
                self.send_header(key, value)
            self.end_headers()
            self.wfile.write(body)

        def _json(self, payload, status: int = 200) -> None:
            self._send(status, json.dumps(payload).encode(), "application/json")

        def _trusted(self) -> bool:
            # Stops other websites reaching this server through the browser
            # (DNS rebinding, cross-site POSTs).
            if self.headers.get("Host") not in allowed_hosts:
                self._json({"error": "Forbidden"}, 403)
                return False
            return True

        def _job_id(self, raw: str) -> str | None:
            if JOB_ID.fullmatch(raw) and (manager.runs_dir / raw / "job.json").exists():
                return raw
            self._json({"error": "No such run."}, 404)
            return None

        def do_GET(self) -> None:
            if not self._trusted():
                return
            parts = [p for p in self.path.split("?")[0].split("/") if p]
            if not parts:
                self._send(200, PAGE.read_bytes(), "text/html; charset=utf-8")
            elif parts == ["api", "jobs"]:
                self._json(manager.list_jobs())
            elif parts == ["api", "presets"]:
                self._json(rg_tools.PRESETS)
            elif len(parts) == 3 and parts[:2] == ["api", "jobs"]:
                if job_id := self._job_id(parts[2]):
                    self._json(manager.detail(job_id))
            elif len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "rg.csv":
                if job_id := self._job_id(parts[2]):
                    path = manager.runs_dir / job_id / "rg.csv"
                    if path.exists():
                        self._send(
                            200,
                            path.read_bytes(),
                            "text/csv",
                            {"Content-Disposition": f'attachment; filename="rg-{job_id}.csv"'},
                        )
                    else:
                        self._json({"error": "No results yet."}, 404)
            else:
                self._json({"error": "Not found"}, 404)

        def do_POST(self) -> None:
            if not self._trusted():
                return
            # A JSON content type can't be sent cross-site without a CORS preflight,
            # which this server never approves.
            if self.headers.get("Content-Type", "").split(";")[0] != "application/json":
                self._json({"error": "Expected JSON."}, 415)
                return
            parts = [p for p in self.path.split("/") if p]
            try:
                length = int(self.headers.get("Content-Length", 0))
                body = json.loads(self.rfile.read(length) or b"{}")
                if not isinstance(body, dict):
                    raise ValueError
            except ValueError:
                self._json({"error": "Invalid JSON."}, 400)
                return

            if parts == ["api", "parse"]:
                try:
                    seq = rg_tools.parse_sequence(str(body.get("sequence", "")))
                except rg_tools.SequenceError as exc:
                    self._json({"error": str(exc)})
                    return
                plans = {p: rg_tools.run_plan(len(seq), p) for p in rg_tools.PRESETS}
                self._json({"sequence": seq, "n_residues": len(seq), "plans": plans})
            elif parts == ["api", "jobs"]:
                try:
                    self._json({"id": manager.submit(validate_job(body))}, 201)
                except BadRequest as exc:
                    self._json({"error": str(exc)}, 400)
            elif len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "cancel":
                if job_id := self._job_id(parts[2]):
                    manager.cancel(job_id)
                    self._json({"ok": True})
            else:
                self._json({"error": "Not found"}, 404)

    return Handler


def serve(port: int) -> tuple[ThreadingHTTPServer, JobManager]:
    manager = JobManager()
    allowed = {f"localhost:{port}", f"127.0.0.1:{port}"}
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(manager, allowed))
    return server, manager
