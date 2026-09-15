import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md
from tools.paths import ensure_runtime_dir
from tools.metadata import write_metadata

def build_sim(sim: Sim):
     components = sim.components

     # Chains sit on an exactly-filled 4 x 4 lattice, one per
     # spacing x spacing cell, and the box is exactly that lattice wide. The
     # periodic images therefore continue the lattice without a seam, so this
     # is an infinite grafted surface and every chain is equally crowded.
     # Placing each chain at the *centre* of its cell (the +0.5) keeps it as far
     # from the box edge as the lattice allows.
     nx, ny = 4, 4
     # spacing = 1/sqrt(0.05 chains/nm^2)
     spacing = 4.472

     ibead = 0
     i = 0
     for comp in components:
          # CALVADOS builds each chain's starting conformation with a lateral
          # offset and extent of its own, so xinit's x/y centre is not (0, 0).
          # Recentre it on the lattice point, otherwise the whole grafting
          # pattern sits offset from the lattice the box is built around.
          xy_centre = 0.5 * (comp.xinit[:, :2].min(axis=0) + comp.xinit[:, :2].max(axis=0))
          xinit_centred = comp.xinit - np.array([xy_centre[0], xy_centre[1], 0.0])

          for idx in range(comp.nmol):
               j = ibead + comp.nbeads

               x0 = (i % nx + 0.5) * spacing
               y0 = (i // nx + 0.5) * spacing

               # Bead 0 is the "Z"-tagged bead. CALVADOS gives "Z" a molecular
               # weight of -2 which the +2 N-terminus patch cancels to exactly
               # 0, and OpenMM holds zero-mass particles completely fixed — so
               # this lattice point is where the chain stays grafted for the
               # whole run, in x, y and z.
               pos = xinit_centred + np.array([x0, y0, 2.0])
               sim.pos[ibead:j] = pos

               ibead = j
               i += 1

     md.Trajectory(sim.pos, sim.top, 0, sim.box, [90,90,90]).save(sim.pdb_cg)


# Job settings for Delft Blue
partition = "gpu-a100"
runtime = "0:30:00"
cpu_per_task = "18"

sim_name = Path(__file__).parent.name

box = [17.888, 17.888, 50.0]
N_save = 7000
N_frames = 1010

# OpenMM runs in nm. The wall keeps the "Z"-tagged end of each ELP anchored
# near the surface at z ~= 0.
z_wall = 1.9

sequences: dict[str, str] = {
     "does-teams-work2": "VPGVGVPGVGVPGVGVPGVGVPGKGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGMGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGKGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGKGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGMGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGVGVPGKGVPGVGVPGVGVPGVGVPGVGVPGVG"
}
sequences["does-teams-work2"] = f"Z{sequences['does-teams-work2'][1:]}"


if __name__ == "__main__":
     path = Path(__file__).parent.resolve()
     cwd = Path(os.getcwd())

     runtime_dir = ensure_runtime_dir(path)

     fasta_file = runtime_dir / "molecules.fasta"

     # Config preparation
     config = Config(
          sysname = sim_name,
          box = box,
          temp = 293.15,
          ionic = 0.19,
          pH = 7.5,
          ext_force = True,
          ext_force_expr = f'step({z_wall}-z)*0.5*({z_wall}-z)^2',
          topol = 'grid',
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
          nmol = 16, # number of molecules
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

     # Settings snapshot next to the trajectory-to-be. run.py rewrites it when
     # the run finishes (status/frame count); `sim metadata` refreshes it any time.
     write_metadata(runtime_dir)
