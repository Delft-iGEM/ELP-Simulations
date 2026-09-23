"""ELP4-30 in free solution — a validation run, not a brush run.

One chain of G(VPGVG)30FC in a cubic box of water-as-continuum, to compare the
model's hydrodynamic radius against the measured value. Nothing here touches
the surface-attached pipeline: this folder only exists to be run and measured.

Three settings differ from every other simulation in this project, and all
three are CALVADOS defaults rather than overrides — which is why free solution
needs no change to tools/:

  topol = 'center'   one molecule at the centre of the box, no grafting lattice
  ext_force = False  no wall; the chain is in bulk, not on a surface
  nmol = 1           a single chain, so nothing it sees is its own neighbour

Residue 0 is NOT tagged "Z". That tag is the surface anchor: CALVADOS gives it
zero mass and OpenMM holds it fixed, which for a free chain would pin one end
and change both the statistics and the molar mass. Here every residue is the
one the sequence says.

Conditions match the experiment being compared against: 20 C, 20 mM salt,
pH 6.5. Note this is NOT the 0.19 M physiological salt the brush runs use; at
20 mM the Debye length is ~2 nm, but ELP4 carries no charged side chains at all
(no K, R, D or E), so only the termini are charged and electrostatics are a
small correction either way.

The box side is set so the chain cannot see its own periodic image: it must
exceed the chain's extent plus the 4 nm electrostatics cutoff.
"""

import os
from pathlib import Path

from calvados.cfg import Config, Components
from calvados.sim import Sim

from tools.paths import ensure_runtime_dir
from tools.metadata import write_metadata

sim_name = Path(__file__).parent.name

sequence = "GVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGFC"          # G(VPGVG)30FC — 153 residues, 12.6 kDa
box_side = 22.0                 # nm, cubic
N_save = 7000
N_frames = 1429
seed = 20260953


def build_sim(sim: Sim):
    """Nothing to do: CALVADOS already placed the single chain at the box centre.

    The brush runs override this hook to lay chains out on a grafting lattice.
    A free chain needs no placement, so this is deliberately empty — and run.py
    calls it either way.
    """
    return


if __name__ == "__main__":
    path = Path(__file__).parent.resolve()
    cwd = Path(os.getcwd())
    runtime_dir = ensure_runtime_dir(path)
    fasta_file = runtime_dir / "molecules.fasta"

    config = Config(
        sysname = sim_name,
        box = [box_side, box_side, box_side],
        temp = 293.15,        # 20 C, the temperature the R_h values were measured at
        ionic = 0.020,        # 20 mM NaCl
        pH = 6.5,
        topol = 'center',     # single chain, centred; no grafting lattice
        ext_force = False,    # no wall: this is bulk solution
        wfreq = N_save,
        steps = N_frames * N_save,
        runtime = 0,
        platform = 'CUDA',
        restart = 'checkpoint',
        frestart = 'restart.chk',
        random_number_seed = seed,
        verbose = True,
    )

    components = Components(
        molecule_type = 'protein',
        nmol = 1,
        restraint = False,
        charge_termini = 'both',    # a real free chain has both termini charged
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
        .replace("{runtime}", "1:00:00").replace("{cpu_per_task}", "18")
    (runtime_dir / "job.sh").write_text(job)

    write_metadata(runtime_dir)
