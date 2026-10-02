"""75th pct of the 24 free-chain RGD arrangements, rerun to 500 ns.

Repeat of simulations/free-rgdarr22-c-heavy — same 352-residue sequence, same box,
same conditions — at 500 ns instead of 200, and on a fresh seed so this is an
independent draw rather than a continuation. The 200 ns sweep put this run
at the 75th percentile (Rg = 5.56 nm over 50-200 ns), but the spread across all 24
arrangements (sd 0.086 nm) was the size of a single run's own error (0.084 nm),
so that ranking was not resolvable. Four runs at 2.5x the length — the two ends
and the two quartiles — test whether the order survives better statistics.

Free solution is CALVADOS's own default, so nothing in tools/ is involved:
topol 'center' places the one molecule at the box centre, ext_force is false so
there is no wall, and residue 0 is not tagged "Z" so nothing is pinned.

Conditions match the grafted runs (293.15 K, 0.19 M, pH 7.5) — the comparison
is the point, so only the geometry may differ.
"""


import os
from pathlib import Path

from calvados.cfg import Config, Components
from calvados.sim import Sim

from tools.paths import ensure_runtime_dir
from tools.metadata import write_metadata

sim_name = Path(__file__).parent.name

sequence = "LQVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGMGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGMGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASA"
box_side = 35.0          # nm, cubic; clears the chain's extent plus the 4 nm cutoff
N_save = 7000
N_frames = 7143
seed = 20260952


def build_sim(sim: Sim):
    """Nothing to place: CALVADOS centres the single chain itself."""
    return


if __name__ == "__main__":
    path = Path(__file__).parent.resolve()
    cwd = Path(os.getcwd())
    runtime_dir = ensure_runtime_dir(path)
    fasta_file = runtime_dir / "molecules.fasta"

    config = Config(
        sysname = sim_name, box = [box_side, box_side, box_side],
        temp = 293.15, ionic = 0.19, pH = 7.5,
        topol = 'center', ext_force = False,
        wfreq = N_save, steps = N_frames * N_save, runtime = 0,
        platform = 'CUDA', restart = 'checkpoint', frestart = 'restart.chk',
        random_number_seed = seed, verbose = True,
    )
    components = Components(
        molecule_type = 'protein', nmol = 1, restraint = False,
        charge_termini = 'both',
        fresidues = str(cwd / "residues_CALVADOS2.csv"),
        ffasta = str(fasta_file),
    )
    components.add(name=sim_name)

    config.write(str(runtime_dir), "config.yaml")
    components.write(str(runtime_dir), "components.yaml")
    fasta_file.write_text(f">{sim_name}\n{sequence}")
    (runtime_dir / "run.py").write_text((cwd / "template/run.py").read_text())
    job = (cwd / "template/job.sh").read_text() \
        .replace("{sim_name}", sim_name).replace("{partition}", "gpu-a100") \
        .replace("{runtime}", "0:50:46").replace("{cpu_per_task}", "18")
    (runtime_dir / "job.sh").write_text(job)
    write_metadata(runtime_dir)
