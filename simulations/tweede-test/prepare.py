import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md

# Parameters
def get_clockwise_position(i):
     x, y = 0, 0
     dx, dy = 0, -1  # Starting direction logic

     for _ in range(i):
          if x == y or (x < 0 and x == -y) or (x > 0 and x == 1 - y):
               # Rotate direction 90 degrees clockwise: (dx, dy) -> (-dy, dx)
               dx, dy = -dy, dx
          x, y = x + dx, y + dy

     return np.array((x, y, 0))

def build_sim(sim: Sim):
     components = sim.components

     ibead = 0
     i = 0
     for comp in components:
          for idx in range(comp.nmol):
               j = ibead + comp.nbeads

               x_c = (sim.box * 0.5)
               x_c[2] = 2

               pos = comp.xinit + x_c + get_clockwise_position(i) * 1.0
               sim.pos[ibead:j] = pos
               ibead = j
               i += 1

     md.Trajectory(sim.pos, sim.top, 0, sim.box, [90,90,90]).save(sim.pdb_cg)


# Job settings for Delft Blue
partition = "gpu-a100"
runtime = "24:30:00"
cpu_per_task = "18"

sim_name = Path(__file__).parent.name

box = [17.32, 17.32, 50.0]
N_save = 5
N_frames = 1000

# OpenMM runs in nm. The wall keeps the "Z"-tagged end of each ELP anchored
# near the surface at z ~= 0.
z_wall = 1.9

sequences: dict[str, str] = {
     "tweede-test": "VPGIGVPGIG"
}
sequences["tweede-test"] = f"Z{sequences['tweede-test'][1:]}"


if __name__ == "__main__":
     path = Path(__file__).parent.resolve()
     cwd = Path(os.getcwd())

     runtime_dir = path / "runtime"
     runtime_dir.mkdir(exist_ok=True)

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
          nmol = 300, # number of molecules
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
