"""even as a single free chain — the solution counterpart of simulations/rgdarr02-even.

Identical 352-residue sequence, identical domain placement. The only difference
is that this one is alone in bulk instead of one of 16 grafted chains. Run as a
pair with the grafted version to separate what the *sequence* does from what the
*layer* does: a domain can be buried by its own chain, which happens in both, or
by its neighbours and the surface, which happens only in the gel.

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

sequence = "LQVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGMGVPGIGVPGIGVPGIGVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGMGVPGIGVPGIGVPGIGVPGIGVPGIG"
box_side = 35.0          # nm, cubic; clears the chain's extent plus the 4 nm cutoff
N_save = 7000
N_frames = 2857
seed = 20260928


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
        .replace("{runtime}", "1:30:00").replace("{cpu_per_task}", "18")
    (runtime_dir / "job.sh").write_text(job)
    write_metadata(runtime_dir)
