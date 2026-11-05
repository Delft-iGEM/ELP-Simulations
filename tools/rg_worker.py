"""Run one CALVADOS simulation for a job folder and write its Rg results.

Usage: python -m tools.rg_worker runs/<job-id>

Reads job.json; writes the CALVADOS files, rg.csv and, last of all, result.json.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import mdtraj as md
import openmm
from calvados.cfg import Components, Config
from calvados.sim import Sim
from yaml import safe_load

from tools import rg as rg_tools

SYSNAME = "elp"
PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESIDUES_CSV = PROJECT_ROOT / "residues_CALVADOS2.csv"


def pick_platform() -> str:
    """Fastest platform that can actually create a context on this machine."""
    for name in ("CUDA", "OpenCL"):
        try:
            system = openmm.System()
            system.addParticle(1.0)
            openmm.Context(
                system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName(name)
            )
            return name
        except Exception as exc:
            print(f"Platform {name} unavailable: {exc}", flush=True)
    return "CPU"


def simulate(job_dir: Path, job: dict, platform: str) -> None:
    plan = job["plan"]
    (job_dir / "molecules.fasta").write_text(f">{SYSNAME}\n{job['sequence']}\n")

    Config(
        sysname=SYSNAME,
        box=[plan["box"]] * 3,
        temp=job["temp_c"] + 273.15,
        ionic=job["ionic"],
        pH=job["ph"],
        topol="center",
        wfreq=plan["wfreq"],
        logfreq=plan["wfreq"],
        steps=plan["steps"],
        runtime=0,
        platform=platform,
        threads=max(1, (os.cpu_count() or 2) // 2),
        restart=None,
        random_number_seed=job["seed"],
        verbose=True,
    ).write(str(job_dir), "config.yaml")

    components = Components(
        molecule_type="protein",
        nmol=1,
        restraint=False,
        charge_termini="both",
        fresidues=str(RESIDUES_CSV),
        ffasta=str(job_dir / "molecules.fasta"),
    )
    components.add(name=SYSNAME)
    components.write(str(job_dir), "components.yaml")

    config = safe_load((job_dir / "config.yaml").read_text())
    components_dict = safe_load((job_dir / "components.yaml").read_text())
    sim = Sim(str(job_dir), config, components_dict)
    sim.build_system()
    sim.simulate()


def analyze(job_dir: Path, job: dict, platform: str) -> None:
    traj = md.load(str(job_dir / f"{SYSNAME}.dcd"), top=str(job_dir / "top.pdb"))
    xyz = rg_tools.unwrap_chain(traj.xyz.astype(float), traj.unitcell_lengths.astype(float))
    rg = rg_tools.radius_of_gyration(xyz, rg_tools.bead_masses(job["sequence"], RESIDUES_CSV))

    summary = rg_tools.summarize(rg)
    frame_ns = job["plan"]["wfreq"] * rg_tools.TIMESTEP_PS / 1000
    skip = summary["frames_discarded"]
    cumulative = rg_tools.cumulative_mean(rg, skip)
    from_start = rg_tools.cumulative_mean(rg)
    with open(job_dir / "rg.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["time_ns", "rg_nm", "rg_cumulative_nm", "rg_cumulative_from_start_nm", "used_in_mean"]
        )
        for i, value in enumerate(rg):
            writer.writerow(
                [
                    f"{(i + 1) * frame_ns:.3f}",
                    f"{value:.4f}",
                    f"{cumulative[i]:.4f}" if i >= skip else "",
                    f"{from_start[i]:.4f}",
                    int(i >= skip),
                ]
            )

    tmp = job_dir / "result.json.tmp"
    tmp.write_text(json.dumps({**summary, "platform": platform, "frame_ns": frame_ns}, indent=2))
    tmp.replace(job_dir / "result.json")


def main() -> None:
    job_dir = Path(sys.argv[1]).resolve()
    job = json.loads((job_dir / "job.json").read_text())
    platform = pick_platform()
    print(f"Using platform {platform}", flush=True)
    simulate(job_dir, job, platform)
    analyze(job_dir, job, platform)
    print("Done", flush=True)


if __name__ == "__main__":
    main()
