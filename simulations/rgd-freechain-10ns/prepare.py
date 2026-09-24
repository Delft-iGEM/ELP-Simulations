"""The project's own RGD-distributed construct as a single free chain — short run.

Same 298-residue sequence as the lysarr sweep (6xHis, 50 VPGIG, 4 VPGKG at
pentapeptide 6/21/35/50, 2 VPGMG at 19/38, 4 RGD after 11/22/33/44), but one
chain in bulk instead of 16 on a surface. A reference for what the molecule
does on its own, and short and finely sampled enough to watch: 1000 frames
over 10 ns is one frame every 10 ps.

Free solution is CALVADOS's own default, so nothing in tools/ is involved:

  topol = 'center'   one molecule at the box centre, no grafting lattice
  ext_force = False  no wall
  nmol = 1           nothing it sees is a neighbour or its own image

Residue 0 is NOT tagged "Z". That anchor has zero mass and is held fixed, which
would pin one end of a free chain.

Conditions follow the *production* runs (0.19 M, pH 7.5, 293.15 K), not the
ELP4 literature-validation runs (0.020 M, pH 6.5) — this is our construct, and
the point is to compare it against our own brush work.
"""

import os
from pathlib import Path

from calvados.cfg import Config, Components
from calvados.sim import Sim

from tools.paths import ensure_runtime_dir
from tools.metadata import write_metadata

sim_name = Path(__file__).parent.name

sequence = "HHHHHHVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGRGDVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGMGVPGIGVPGKGVPGIGRGDVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGRGDVPGIGVPGKGVPGIGVPGIGVPGMGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGRGDVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIG"
box_side = 30.0          # nm, cubic — clears the chain plus the 4 nm cutoff
N_save = 1000          # steps per frame: 1000 x 0.01 ps = 10 ps
N_frames = 1000         # -> 10 ns in 1000 frames
seed = 20260924


def build_sim(sim: Sim):
    """Nothing to place: CALVADOS centres the single chain itself."""
    return


if __name__ == "__main__":
    path = Path(__file__).parent.resolve()
    cwd = Path(os.getcwd())
    runtime_dir = ensure_runtime_dir(path)
    fasta_file = runtime_dir / "molecules.fasta"

    config = Config(
        sysname = sim_name,
        box = [box_side, box_side, box_side],
        temp = 293.15,
        ionic = 0.19,          # as the brush runs, not the ELP4 validation runs
        pH = 7.5,
        topol = 'center',
        ext_force = False,
        wfreq = N_save,
        steps = N_frames * N_save,
        runtime = 0,
        platform = 'CPU',
        threads = 4,
        restart = 'checkpoint',
        frestart = 'restart.chk',
        random_number_seed = seed,
        verbose = True,
    )

    components = Components(
        molecule_type = 'protein',
        nmol = 1,
        restraint = False,
        charge_termini = 'both',
        fresidues = str(cwd / "residues_CALVADOS2.csv"),
        ffasta = str(fasta_file),
    )
    components.add(name=sim_name)

    config.write(str(runtime_dir), "config.yaml")
    components.write(str(runtime_dir), "components.yaml")
    fasta_file.write_text(f">{sim_name}\n{sequence}")
    (runtime_dir / "run.py").write_text((cwd / "template/run.py").read_text())
    write_metadata(runtime_dir)
