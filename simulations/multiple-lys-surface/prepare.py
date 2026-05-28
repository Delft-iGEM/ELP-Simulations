import os

import numpy as np
import mdtraj as md

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path

from tools.elibpy import build_sequence_with_features

# Parameters
# OpenMM runs in nm where as PDB is in Angtrums so divide it by 10
z_wall = 230 / 10

def build_sim(sim: Sim):
    # Translate the ENTIRE chain so bead 0 (Head-K / first X) lands at z_wall.
    #
    # Why not just set x_bead_indices[:,2] = z_wall?
    # build_linear() grows the chain in +z centered at 0, so for 303 beads
    # bead 0 starts at z ≈ -57 nm → after topol='center' it's at z ≈ -32 nm.
    # Jumping only the X beads to z=23 while their neighbours stay at -32 nm
    # creates ~55 nm bond stretches → forces blow up → NaN.
    # Shifting every bead by the same Δz keeps all bonds intact.
    pos_array = np.array(sim.pos)          # (n_atoms, 3) in nm

    z_shift = z_wall - pos_array[0, 2]    # bead 0 = Head-K anchor
    pos_array[:, 2] += z_shift
    sim.pos = list(pos_array)

    traj = md.Trajectory(
        pos_array.reshape(1, -1, 3), sim.top,
        unitcell_lengths=[float(sim.box[0]), float(sim.box[1]), float(sim.box[2])],
        unitcell_angles=[90.0, 90.0, 90.0],
    )
    traj.save_pdb(sim.pdb_cg)

# Job settings for Delft Blue
partition = "gpu-a100-small"
runtime = "3:30:00"
cpu_per_task = "2"

# or
# partition = "gpu-a100"
# runtime = "24:30:00"
# cpu_per_task = "18"

sim_name = Path(__file__).parent.name

L = 50
N_save = 7000
N_frames = 1010

sequences: dict[str, str] = {
     "K-LDAS-RGD-SA-I2-K1-I12-VP": build_sequence_with_features([
          ("K", "Head"),
          ("LDAS", "Buffer"),
          ("TVYAVTGRGDSPASSAA", "RGD"),
          ("SA", "Buffer"),
          ("I2", "I"),
          ("K1", "K"),
          ("I12", "I"),
          ("VP", "Buffer")
     ])["seq"] * 3
}

sequences["K-LDAS-RGD-SA-I2-K1-I12-VP"] = sequences["K-LDAS-RGD-SA-I2-K1-I12-VP"].replace("K", "X")


if __name__ == "__main__":
     path = Path(__file__).parent.resolve()
     cwd = Path(os.getcwd())

     runtime_dir = path / "runtime"
     runtime_dir.mkdir(exist_ok=True)

     fasta_file = runtime_dir / "molecules.fasta"

     # Config preparation
     config = Config(
          sysname = sim_name,
          box = [L, L, L],
          temp = 293.15,
          ionic = 0.19,
          pH = 7.5,
          ext_force = True,
          ext_force_expr = f'step({z_wall}-z)*0.5*({z_wall}-z)^2',
          topol = 'center',
          wfreq = N_save,
          steps = N_frames*N_save,
          runtime = 0,
          platform = 'CUDA',
          restart = 'checkpoint',
          frestart = 'restart.chk',
          verbose = True
     )

     components = Components(
          # Defaults
          molecule_type = 'protein',
          nmol = 1, # number of molecules
          restraint = False,
          charge_termini = 'both',
          fresidues = str(cwd / "residues_CALVADOS2.csv"),
          ffasta = str(fasta_file)
     )

     for k in sequences.keys():
          components.add(name=k)

     # Prepare runtime - You probably don't need to edit below this

     # Write config and components
     config.write(str(runtime_dir), "config.yaml")
     components.write(str(runtime_dir), "components.yaml")

     # Create FASTA
     fasta_content = "\n".join([f">{k}\n{v}" for k, v in sequences.items()])
     
     fasta_file.write_text(fasta_content)
     
     # Move job.sh and run.py
     run_file = (cwd / "template/run.py").read_text()
     (runtime_dir / "run.py").write_text(run_file)

     job_file = (cwd / "template/job.sh").read_text()
     
     job_file = job_file \
          .replace("{sim_name}", sim_name) \
          .replace("{partition}", partition) \
          .replace("{runtime}", runtime) \
          .replace("{cpu_per_task}", cpu_per_task)

     (runtime_dir / "job.sh").write_text(job_file)