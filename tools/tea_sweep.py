"""Sweep driver for CALVADOS-TEA: sequences x temperatures, one run per cell.

    .venv/bin/python tools/tea_sweep.py init          # write the workspace
    .venv/bin/python tools/tea_sweep.py submit        # sbatch everything not yet queued
    .venv/bin/python tools/tea_sweep.py status

Layout is signac-compatible: ``workspace/<job_id>/`` with a ``signac_statepoint.json``
per job, where job_id is the md5 of the canonical state point, so ``signac.get_project()``
picks the workspace up without this module being installed. signac itself is not
imported — it is not a dependency of this project and the driver works without it.

Each job directory carries a manifest with everything needed to reproduce the
cell: gamma, T0, the full lambda vector as simulated, salt, box, seed, CALVADOS
version and the git hash of this tree.

The free single chains here are NOT the surface pipeline: no wall, no anchor,
one molecule at the box centre, exactly like simulations/free-rgdarr*.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from tools.new_simulation import _project_root
from tools.paths import data_root
from tools.tea import DEFAULT_GAMMA, T0, lambda_table, write_residues_csv

# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------

SEQUENCES = {               # (VPGXG)84 -> 420 residues
    "V84": "VPGVG" * 84,
    "I84": "VPGIG" * 84,
    "F84": "VPGFG" * 84,
    "E84": "VPGEG" * 84,
}
TEMPERATURES = [300.0, 310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0]
GAMMA = DEFAULT_GAMMA

# Authors' single-chain settings, except where their LAMMPS setup does not carry
# over to CALVADOS 2. Deviations from stock are recorded in each manifest.
BOX_NM = 80.0
TIMESTEP_PS = 0.01          # 10 fs

# Stock CALVADOS. The paper's 0.5/ps is from their LAMMPS runs; equilibrium Rg
# does not depend on friction, so 0.5 would only slow configurational sampling.
FRICTION_PER_PS = 0.01

# Stock CALVADOS 2, NOT the paper's 4 nm. CALVADOS 2 is the tuned-interaction-
# range model: its lambda(T0) values — the anchor the whole of TEA is expressed
# relative to — were parameterised at a 2 nm LJ cutoff. Forcing 4 nm would throw
# that parameterisation away while still using the lambdas it produced. The
# paper's 4 nm goes with HPS-Urry lambdas, which this pipeline does not use.
CUTOFF_LJ_NM = 2.0
CUTOFF_YU_NM = 4.0          # electrostatics, stock

# Deliberately the paper's value, not this project's stock 0.19 M, so the sweep
# is comparable with their published numbers. Flagged in every manifest.
IONIC_M = 0.1
STOCK_IONIC_M = 0.19
PH = 7.5
PRODUCTION_US = 4.0         # >200 residues
DISCARD_NS = 50.0
SAVE_PS = 100.0             # frame every 100 ps -> 40 000 frames over 4 us

# Measured, not guessed: tools/tea_sweep.py benchmark gave 974.4 ns/h (27,068
# steps/s) for this exact system — 420 beads, 80 nm box, 2/4 nm cutoffs, A100.
# 4 us is therefore ~4.1 h; 6 h is ~1.5x headroom and well under DelftBlue's
# 24 h per-job cap. Re-run the benchmark if the box, cutoffs or chain change.
MEASURED_NS_PER_HOUR = 974.4
WALLTIME = "6:00:00"
PARTITION = "gpu-a100"
CPUS = "18"
SEED0 = 20261002

WORKSPACE = "tea-sweep"

# One short run at the cheapest point of the sweep, to measure ns/hour on the
# real system before committing 32 jobs. 420 beads is what the sweep runs; a
# steps/s figure from a different system size or from a slab does not transfer.
BENCHMARK_NAME = "tea-benchmark"
BENCHMARK_NS = 10.0
BENCHMARK_WALLTIME = "0:30:00"


def parse_hours(walltime: str) -> float:
    h, m, sec = (int(x) for x in walltime.split(":"))
    return h + m / 60.0 + sec / 3600.0


def _steps_and_save() -> tuple[int, int]:
    n_save = int(round(SAVE_PS / TIMESTEP_PS))
    n_frames = int(round(PRODUCTION_US * 1e6 / SAVE_PS))
    return n_save, n_frames


def state_point(sequence_name: str, temperature: float) -> dict:
    return {"sequence": sequence_name, "temperature_K": temperature,
            "gamma": GAMMA, "model": "CALVADOS2-TEA", "ensemble": "single-chain"}


def job_id(sp: dict) -> str:
    canonical = json.dumps(sp, sort_keys=True, separators=(",", ":"))
    return hashlib.md5(canonical.encode()).hexdigest()


def _git_hash(root: Path) -> str:
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True)
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                               capture_output=True, text=True).stdout.strip()
        return out.stdout.strip() + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def _calvados_version() -> str:
    try:
        import calvados
        return getattr(calvados, "__version__", "unknown")
    except Exception:
        return "unknown"


PREPARE = '''"""(VPGXG)84 free single chain at one temperature, CALVADOS-TEA.

Generated by tools/tea_sweep.py — edit the driver, not this file.
"""

import os
from pathlib import Path

from calvados.cfg import Config, Components
from calvados.sim import Sim

from tools.paths import ensure_runtime_dir
from tools.metadata import write_metadata
from tools.tea import report as tea_report, write_residues_csv

sim_name = Path(__file__).parent.name

sequence = "{sequence}"
temperature = {temperature!r}
gamma = {gamma!r}
box_side = {box!r}
N_save = {n_save}
N_frames = {n_frames}
seed = {seed}


def build_sim(sim: Sim):
    """CALVADOS centres the single chain itself."""
    return


if __name__ == "__main__":
    path = Path(__file__).parent.resolve()
    cwd = Path(os.getcwd())
    runtime_dir = ensure_runtime_dir(path)
    fasta_file = runtime_dir / "molecules.fasta"

    # TEA: lambda at this run's temperature, written here and used only here.
    # The project's residues_CALVADOS2.csv is read and never modified.
    fresidues = runtime_dir / "residues_TEA.csv"
    lambdas = write_residues_csv(fresidues, temperature, gamma,
                                 source=cwd / "residues_CALVADOS2.csv")
    (runtime_dir / "tea.yaml").write_text(
        "# CALVADOS-TEA, as this run used it. Read back by tools/metadata.py.\\n"
        f"temperature_K: {{temperature}}\\n"
        f"gamma: {{gamma}}\\n"
        "T0_K: 300.0\\n"
        "source_table: residues_CALVADOS2.csv\\n"
        "lambdas:\\n"
        + "".join(f"  {{k}}: {{v!r}}\\n" for k, v in sorted(lambdas.items())))
    print(tea_report(temperature, gamma))

    config = Config(
        sysname = sim_name, box = [box_side, box_side, box_side],
        temp = temperature, ionic = {ionic!r}, pH = {ph!r},
        topol = 'center', ext_force = False,
        cutoff_lj = {cutoff_lj!r}, cutoff_yu = {cutoff_yu!r},
        friction_coeff = {friction!r},
        wfreq = N_save, steps = N_frames * N_save, runtime = 0,
        platform = 'CUDA', restart = 'checkpoint', frestart = 'restart.chk',
        random_number_seed = seed, verbose = True,
    )
    components = Components(
        molecule_type = 'protein', nmol = 1, restraint = False,
        charge_termini = 'both',
        fresidues = str(fresidues),
        ffasta = str(fasta_file),
    )
    components.add(name=sim_name)

    config.write(str(runtime_dir), "config.yaml")
    components.write(str(runtime_dir), "components.yaml")
    fasta_file.write_text(f">{{sim_name}}\\n{{sequence}}")
    (runtime_dir / "run.py").write_text((cwd / "template/run.py").read_text())
    job = (cwd / "template/job.sh").read_text() \\
        .replace("{{sim_name}}", sim_name).replace("{{partition}}", "{partition}") \\
        .replace("{{runtime}}", "{walltime}").replace("{{cpu_per_task}}", "{cpus}")
    (runtime_dir / "job.sh").write_text(job)
    write_metadata(runtime_dir)
'''


def init(root: Path, dry_run: bool = False) -> list[tuple[str, Path]]:
    """Write one simulation folder + signac state point per (sequence, T)."""
    n_save, n_frames = _steps_and_save()
    workspace = root / WORKSPACE / "workspace"
    made = []
    for s_index, (seq_name, sequence) in enumerate(sorted(SEQUENCES.items())):
        for t_index, temperature in enumerate(TEMPERATURES):
            sp = state_point(seq_name, temperature)
            jid = job_id(sp)
            name = f"tea-{seq_name.lower()}-{int(temperature)}k"
            sim_dir = root / "simulations" / name
            job_dir = workspace / jid
            made.append((name, sim_dir))
            if dry_run:
                continue

            sim_dir.mkdir(parents=True, exist_ok=True)
            (sim_dir / "prepare.py").write_text(PREPARE.format(
                sequence=sequence, temperature=temperature, gamma=GAMMA,
                box=BOX_NM, n_save=n_save, n_frames=n_frames,
                seed=SEED0 + s_index * 100 + t_index,
                ionic=IONIC_M, ph=PH, cutoff_lj=CUTOFF_LJ_NM, cutoff_yu=CUTOFF_YU_NM,
                friction=FRICTION_PER_PS, partition=PARTITION,
                walltime=WALLTIME, cpus=CPUS))

            job_dir.mkdir(parents=True, exist_ok=True)
            (job_dir / "signac_statepoint.json").write_text(json.dumps(sp, indent=2))
            (job_dir / "manifest.json").write_text(json.dumps({
                "state_point": sp,
                "simulation": name,
                "sim_dir": str(sim_dir),
                "runtime_dir": str(data_root() / name),
                "sequence": sequence,
                "n_residues": len(sequence),
                "tea": {"gamma": GAMMA, "T0_K": T0,
                        "lambdas": lambda_table(temperature, GAMMA)},
                "simulation_settings": {
                    "box_nm": [BOX_NM] * 3, "ionic_M": IONIC_M, "pH": PH,
                    "timestep_ps": TIMESTEP_PS, "friction_per_ps": FRICTION_PER_PS,
                    "cutoff_lj_nm": CUTOFF_LJ_NM, "cutoff_yu_nm": CUTOFF_YU_NM,
                    "steps": n_save * n_frames, "save_every_steps": n_save,
                    "production_us": PRODUCTION_US, "discard_ns": DISCARD_NS,
                    "seed": SEED0 + s_index * 100 + t_index,
                },
                "deviations_from_stock": {
                    "ionic_M": {"used": IONIC_M, "stock": STOCK_IONIC_M,
                                "why": "the paper's value, for comparability with their "
                                       "published single-chain dimensions"},
                },
                "calvados_version": _calvados_version(),
                "git_hash": _git_hash(root),
                "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            }, indent=2))

    if not dry_run:
        (root / WORKSPACE / "signac.rc").write_text(f"project = {WORKSPACE}\n")
        (root / WORKSPACE / ".signac" / "config").parent.mkdir(exist_ok=True)
        (root / WORKSPACE / ".signac" / "config").write_text(f'project = "{WORKSPACE}"\n')
    return made


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command",
                        choices=["init", "submit", "status", "plan", "benchmark"])
    parser.add_argument("--yes", action="store_true",
                        help="submit without the confirmation prompt")
    args = parser.parse_args(argv)
    root = _project_root()

    if args.command == "benchmark":
        n_save, _ = _steps_and_save()
        n_frames = int(round(BENCHMARK_NS * 1000.0 / SAVE_PS))
        sim_dir = root / "simulations" / BENCHMARK_NAME
        sim_dir.mkdir(parents=True, exist_ok=True)
        (sim_dir / "prepare.py").write_text(PREPARE.format(
            sequence=SEQUENCES["V84"], temperature=300.0, gamma=GAMMA,
            box=BOX_NM, n_save=n_save, n_frames=n_frames, seed=SEED0,
            ionic=IONIC_M, ph=PH, cutoff_lj=CUTOFF_LJ_NM, cutoff_yu=CUTOFF_YU_NM,
            friction=FRICTION_PER_PS, partition=PARTITION,
            walltime=BENCHMARK_WALLTIME, cpus=CPUS))
        print(f"wrote simulations/{BENCHMARK_NAME}/ — (VPGVG)84, {len(SEQUENCES['V84'])} beads, "
              f"300 K, {BENCHMARK_NS:g} ns ({n_save * n_frames:,} steps), "
              f"same box/cutoffs/friction as the sweep")
        print(f"  sim prepare {BENCHMARK_NAME} && sim submit {BENCHMARK_NAME}")
        print(f"  then read ns_per_hour / steps_per_second from its metadata.csv")
        return 0

    if args.command in ("plan", "init"):
        made = init(root, dry_run=(args.command == "plan"))
        n_save, n_frames = _steps_and_save()
        print(f"{len(made)} runs: {len(SEQUENCES)} sequences x {len(TEMPERATURES)} temperatures")
        print(f"  {n_save * n_frames:,} steps = {PRODUCTION_US:g} us, {n_frames:,} frames "
              f"at {SAVE_PS:g} ps, discard first {DISCARD_NS:g} ns")
        print(f"  box {BOX_NM:g} nm, cutoffs {CUTOFF_LJ_NM:g}/{CUTOFF_YU_NM:g} nm "
              f"(stock CALVADOS 2), friction {FRICTION_PER_PS:g}/ps (stock), "
              f"gamma {GAMMA:g}, walltime {WALLTIME} each")
        print(f"  ionic {IONIC_M:g} M — deliberately the paper's value, not this project's "
              f"stock {STOCK_IONIC_M:g} M; recorded in every manifest")
        hours = PRODUCTION_US * 1000.0 / MEASURED_NS_PER_HOUR
        print(f"  measured {MEASURED_NS_PER_HOUR:g} ns/h -> ~{hours:.1f} h per run, "
              f"{len(made) * hours:.0f} GPU-hours for the sweep "
              f"(walltime {WALLTIME} = {parse_hours(WALLTIME) / hours:.1f}x headroom)")
        if args.command == "plan":
            print("\n(plan only — nothing written; run `init` to write them)")
        else:
            print(f"\nwrote simulations/ folders and {WORKSPACE}/workspace/<id>/")
            print(f"next: prepare each with `sim prepare <name>`, then "
                  f"`{Path(sys.argv[0]).name} submit`")
        return 0

    if args.command == "status":
        for name, sim_dir in init(root, dry_run=True):
            runtime = sim_dir / "runtime"
            meta = runtime / "metadata.csv"
            state = "not prepared"
            if meta.is_file():
                rows = dict(line.split(",")[:2] for line in meta.read_text().splitlines()[1:]
                            if "," in line)
                state = rows.get("status", "prepared")
            print(f"  {name:<24} {state}")
        return 0

    # submit
    runs = init(root, dry_run=True)
    pending = [(n, d) for n, d in runs if (d / "runtime" / "job.sh").is_file()]
    missing = [n for n, d in runs if not (d / "runtime" / "job.sh").is_file()]
    if missing:
        print(f"{len(missing)} run(s) not prepared yet — run `sim prepare <name>` first:")
        for name in missing[:5]:
            print(f"  {name}")
        return 1
    print(f"about to sbatch {len(pending)} jobs at {WALLTIME} each "
          f"({len(pending) * 24} GPU-hours requested)")
    if not args.yes:
        print("re-run with --yes to actually submit")
        return 0
    for name, sim_dir in pending:
        subprocess.run(["sbatch", str(sim_dir / "runtime" / "job.sh")], cwd=root, check=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
