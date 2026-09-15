import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md
from tools.paths import ensure_runtime_dir
from tools.metadata import write_lattice, write_metadata
from tools.crosslink import CrosslinkSettings, write_settings

# Grafting lattice: 4 x 4 points, one chain per spacing x spacing
# cell, with an empty `margin` rim between the outermost points and the box
# wall on every side (box = nx*spacing + 2*margin). The rim keeps chains away
# from the periodic boundary, so they rarely lean across it and get wrapped to
# the far side of the trajectory; the cost is that neighbours across the
# boundary sit spacing + 2*margin apart, so rim chains are a little less
# crowded than the ones in the middle. margin = 0 restores the seamless tiling.
nx, ny = 4, 4
# spacing = 1/sqrt(0.05 chains/nm^2)
spacing = 4.472
# margin = 0: the lattice fills the box exactly, so the periodic images continue it seamlessly and every chain is equally crowded
margin = 0.0

# Reactive crosslinking. None = off, and then run.py takes the stock CALVADOS
# path unchanged. A dict turns it on — see tools/crosslink.py for what the
# numbers mean and, just as importantly, what they don't.
# lysine pairs within 1 nm (10.0 A) bond during the run, valence 1, p=1.0
crosslink = {'distance': 1, 'valence': 1, 'prob': 1.0, 'check_every': 1000, 'k': 2000.0, 'r0': 0.6, 'ramp_steps': 500, 'selection': None, 'start_step': 0, 'seed': None}


def build_sim(sim: Sim):
     components = sim.components

     # Each chain goes at the *centre* of its own cell (the +0.5), offset by the
     # margin — exactly what tools.new_simulation.lattice_positions previews.

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

               x0 = margin + (i % nx + 0.5) * spacing
               y0 = margin + (i // nx + 0.5) * spacing

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
runtime = "0:15:00"
cpu_per_task = "18"

sim_name = Path(__file__).parent.name

box = [17.888, 17.888, 50.0]
N_save = 7000
N_frames = 1010

# OpenMM runs in nm. The wall keeps the "Z"-tagged end of each ELP anchored
# near the surface at z ~= 0.
z_wall = 1.9

sequences: dict[str, str] = {
     "does-link-work1": "VPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGRGDVPGIGVPGIGVPGMGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGRGDVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIG"
}
sequences["does-link-work1"] = f"Z{sequences['does-link-work1'][1:]}"


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

     # The box alone doesn't say where the lattice inside it sits (spacing and
     # margin are two unknowns in one number), so record them next to the
     # trajectory for metadata.csv to pick up.
     write_lattice(runtime_dir, nx=nx, ny=ny, spacing=spacing, margin=margin)

     # run.py reads this file to decide whether the run reacts. Written (or
     # removed) every time, so a regenerated non-reactive run can never inherit
     # a stale crosslink.yaml from a previous attempt.
     write_settings(runtime_dir,
                    CrosslinkSettings.from_dict(crosslink) if crosslink else None)

     # Settings snapshot next to the trajectory-to-be. run.py rewrites it when
     # the run finishes (status/frame count); `sim metadata` refreshes it any time.
     write_metadata(runtime_dir)
