"""Local web UI: paste an ELP sequence, simulate it, read off the Rg.

Started with `sim ui`. Listens on localhost only. Jobs run one at a time in a
separate process (tools.rg_worker), each in its own folder under runs/.
"""

from __future__ import annotations

import base64
import binascii
import csv
import json
import re
import secrets
import subprocess
import sys
import threading
import unicodedata
from collections import deque
from datetime import datetime
from functools import cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np

from tools import rg as rg_tools

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = PROJECT_ROOT / "runs"
PAGE = Path(__file__).with_name("ui.html")
RACE_SCRIPT = Path(__file__).with_name("race.js")
# Pictures behind the race chart's Rg regions, served at /images/<name>. A fixed list, so no other file is reachable.
REGION_IMAGES = {name: PROJECT_ROOT / "Images" / name for name in ("Good_bone.png", "neuron.png", "skin.png")}
TEAMS_CSV = PROJECT_ROOT / "teams-2026.csv"  # iGEM's 2026 team export
JOB_ID = re.compile(r"\d{8}-\d{6}(-[a-z0-9-]+)?")
# Every run uses the same conditions, so runs compare on sequence alone.
CONDITIONS = {"temp_c": 20.0, "ionic": 0.13, "ph": 7.5, "preset": "fast"}
TARGET_RG_LIMITS = (0.5, 30.0)  # nm
ERROR_TAIL_LINES = 25
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
PNG_DATA_URL = "data:image/png;base64,"
MAX_IMAGE_BYTES = 8 * 1024 * 1024
MAX_BODY_BYTES = 2 * MAX_IMAGE_BYTES  # base64 inflates the image by a third


class BadRequest(ValueError):
    pass


def team_key(name: str) -> str:
    """Spelling-insensitive key: 'Ashesi Ghana', 'ashesi-ghana' and 'AshesiGhana' all match.

    The page's teamKey() must do the same.
    """
    ascii_name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]", "", ascii_name.lower())


@cache
def load_teams() -> list[dict]:
    """Accepted 2026 teams, sorted by name. Empty if the CSV is missing."""
    try:
        with open(TEAMS_CSV, newline="", encoding="utf-8-sig") as f:
            rows = [r for r in csv.DictReader(f) if r.get("status") == "accepted"]
    except OSError:
        return []
    return sorted(({"name": r["name"], "city": r["city"]} for r in rows), key=lambda t: t["name"].lower())


def canonical_team(raw: object) -> str | None:
    text = str(raw or "").strip()
    if not text:
        return None
    for team in load_teams():
        if team_key(team["name"]) == team_key(text):
            return team["name"]
    raise BadRequest(f"'{text[:40]}' isn't a 2026 iGEM team. Pick one from the list.")


def validate_job(body: dict) -> dict:
    try:
        sequence = rg_tools.parse_sequence(str(body.get("sequence", "")))
    except rg_tools.SequenceError as exc:
        raise BadRequest(str(exc)) from exc

    low, high = TARGET_RG_LIMITS
    try:
        target = float(body["target_rg"])
    except (KeyError, TypeError, ValueError):
        raise BadRequest("Enter a target Rg in nm.") from None
    if not low <= target <= high:
        raise BadRequest(f"Target Rg must be between {low:g} and {high:g} nm.")

    return {
        "name": str(body.get("name") or "").strip()[:80],
        "team": canonical_team(body.get("team")),
        "input": str(body["sequence"]).strip(),
        "sequence": sequence,
        "plan": rg_tools.run_plan(len(sequence), CONDITIONS["preset"]),
        "seed": secrets.randbelow(2**31 - 1) + 1,
        "target_rg": target,
        **CONDITIONS,
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
            "team": job.get("team"),  # absent on runs made before teams
            "input": job["input"],
            "n_residues": len(job["sequence"]),
            "temp_c": job["temp_c"],
            "ionic": job["ionic"],
            "ph": job["ph"],
            "preset": job["preset"],
            "simulated_ns": job["plan"]["simulated_ns"],
            "has_image": (job_dir / "image.png").exists(),
            "target_rg": job.get("target_rg"),  # absent on runs made before targets
            "blocks": rg_tools.analyze_blocks(job["input"]),
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

    def set_image(self, job_id: str, data_url: str | None) -> None:
        """Store a run's picture, or remove it when data_url is None.

        The page converts every picture to PNG first, so only PNG is accepted here.
        """
        path = self.runs_dir / job_id / "image.png"
        if data_url is None:
            path.unlink(missing_ok=True)
            return
        if not isinstance(data_url, str) or not data_url.startswith(PNG_DATA_URL):
            raise BadRequest("Expected a PNG picture.")
        try:
            image = base64.b64decode(data_url[len(PNG_DATA_URL) :], validate=True)
        except binascii.Error:
            raise BadRequest("That picture couldn't be read.") from None
        if not image.startswith(PNG_SIGNATURE):
            raise BadRequest("Expected a PNG picture.")
        if len(image) > MAX_IMAGE_BYTES:
            raise BadRequest("That picture is too large.")
        tmp = path.with_suffix(f".{threading.get_ident()}.tmp")
        tmp.write_bytes(image)
        tmp.replace(path)

    def coords(self, job_id: str) -> bytes:
        """COM-centred chain coordinates as little-endian float32 (frames x beads x 3, nm).

        Cached beside the trajectory after the first request.
        """
        job_dir = self.runs_dir / job_id
        cache = job_dir / "coords.bin"
        if not cache.exists():
            job = json.loads((job_dir / "job.json").read_text())
            xyz = rg_tools.load_chain(job_dir / "elp.dcd", job_dir / "top.pdb")
            masses = rg_tools.bead_masses(job["sequence"], rg_tools.RESIDUES_CSV)
            xyz -= np.average(xyz, axis=1, weights=masses)[:, None, :]
            tmp = job_dir / f"coords.{threading.get_ident()}.tmp"
            tmp.write_bytes(xyz.astype("<f4").tobytes())
            tmp.replace(cache)
        return cache.read_bytes()

    def detail(self, job_id: str) -> dict:
        out = self.summary(job_id)
        job_dir = self.runs_dir / job_id
        job = json.loads((job_dir / "job.json").read_text())
        out["sequence"] = job["sequence"]
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
            elif parts == ["race.js"]:
                self._send(200, RACE_SCRIPT.read_bytes(), "text/javascript; charset=utf-8")
            elif len(parts) == 2 and parts[0] == "images" and parts[1] in REGION_IMAGES:
                path = REGION_IMAGES[parts[1]]
                if path.exists():
                    self._send(200, path.read_bytes(), "image/png", {"X-Content-Type-Options": "nosniff"})
                else:
                    self._json({"error": "No such picture."}, 404)
            elif parts == ["api", "jobs"]:
                self._json(manager.list_jobs())
            elif parts == ["api", "teams"]:
                self._json(load_teams())
            elif parts == ["api", "settings"]:
                preset = CONDITIONS["preset"]
                self._json({
                    **CONDITIONS,
                    "preset_label": rg_tools.PRESETS[preset]["label"],
                    "target_rg_limits": TARGET_RG_LIMITS,
                    "elp_guests": rg_tools.ELP_GUESTS,
                    "elp_repeats": rg_tools.ELP_REPEATS,
                    "named_blocks": rg_tools.NAMED_BLOCKS,
                    "block_aliases": rg_tools.BLOCK_ALIASES,
                })
            elif len(parts) == 3 and parts[:2] == ["api", "jobs"]:
                if job_id := self._job_id(parts[2]):
                    self._json(manager.detail(job_id))
            elif len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "image":
                if job_id := self._job_id(parts[2]):
                    path = manager.runs_dir / job_id / "image.png"
                    if path.exists():
                        self._send(
                            200,
                            path.read_bytes(),
                            "image/png",
                            {"X-Content-Type-Options": "nosniff"},
                        )
                    else:
                        self._json({"error": "No picture for this run."}, 404)
            elif len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "coords":
                if job_id := self._job_id(parts[2]):
                    try:
                        self._send(200, manager.coords(job_id), "application/octet-stream")
                    except OSError:
                        self._json({"error": "No trajectory for this run."}, 404)
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
                if length > MAX_BODY_BYTES:
                    self._json({"error": "That picture is too large."}, 413)
                    return
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
                plan = rg_tools.run_plan(len(seq), CONDITIONS["preset"])
                blocks = rg_tools.analyze_blocks(str(body.get("sequence", "")))
                self._json({"sequence": seq, "n_residues": len(seq), "plan": plan, "blocks": blocks})
            elif parts == ["api", "jobs"]:
                try:
                    self._json({"id": manager.submit(validate_job(body))}, 201)
                except BadRequest as exc:
                    self._json({"error": str(exc)}, 400)
            elif len(parts) == 4 and parts[:2] == ["api", "jobs"] and parts[3] == "image":
                if job_id := self._job_id(parts[2]):
                    try:
                        manager.set_image(job_id, body.get("data"))
                        self._json({"ok": True})
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
