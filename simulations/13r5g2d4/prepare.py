import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path
import numpy as np
import mdtraj as md
from tools.paths import ensure_runtime_dir
from tools.metadata import write_lattice, write_metadata
from tools.crosslink import CrosslinkSettings, write_settings
from tools.surface import SurfaceSettings, place_chains, write_settings as write_surface_settings

# Grafting lattice: 4 x 4 points, one chain per spacing x spacing
# cell, with an empty `margin` rim between the outermost points and the box
# wall on every side (box = nx*spacing + 2*margin). The rim keeps chains away
# from the periodic boundary, so they rarely lean across it and get wrapped to
# the far side of the trajectory; the cost is that neighbours across the
# boundary sit spacing + 2*margin apart, so rim chains are a little less
# crowded than the ones in the middle. margin = 0 restores the seamless tiling.
nx, ny = 4, 4
# spacing = 1/sqrt(0.04 chains/nm^2)
spacing = 5.0
# margin = 0: the lattice fills the box exactly, so the periodic images continue it seamlessly and every chain is equally crowded
margin = 0.0

# Reactive crosslinking. None = off, and then run.py takes the stock CALVADOS
# path unchanged. A dict turns it on — see tools/crosslink.py for what the
# numbers mean and, just as importantly, what they don't.
# lysine pairs within 0.8 nm (8.0 A) bond during the run, valence 1, p=1.0
crosslink = {'distance': 0.8, 'valence': 1, 'prob': 1.0, 'check_every': 1000, 'k': 2000.0, 'r0': 0.6, 'ramp_steps': 500, 'selection': None, 'start_step': 1000000, 'seed': None, 'min_span': 3}

# How the chains are attached to the surface. "brush" (the default, and the
# pipeline as it always was): residue 0 of every chain is tagged "Z" and pinned
# to its lattice point. "free" / "preattached": residue 0 is an ordinary
# residue and each chain hangs from lysines bonded to the surface instead —
# the dict below holds those settings, and tools/surface.py says what they
# mean and what the results can and cannot be trusted for.
# preattached: lysines within 0.8 nm (8.0 A) of the tether plane bond to it; bound lysines sit 0.6 nm above the wall; cap 0.5 of all lysines; 0.3 bound at t=0
mode = "preattached"
surface = {'mode': 'preattached', 'distance': 0.8, 'prob': 1.0, 'max_fraction': 0.5, 'preattached_fraction': 0.3, 'dynamic': None, 'k': 2000.0, 'tether': 0.6, 'ramp_steps': 500, 'check_every': 1000, 'start_step': 1000000, 'attraction': 0.0, 'attraction_width': 0.5, 'seed': 2023927060, 'z_wall': 1.9}


def build_sim(sim: Sim):
     components = sim.components

     if surface is not None:
          # free / preattached: every chain stands on its pinned lysine(s).
          # tools.surface builds the starting conformations on the same lattice
          # cells as below, verifies them, and leaves the pin plan on `sim` for
          # run.py's reaction loop.
          place_chains(sim, SurfaceSettings.from_dict(surface),
                       nx=nx, ny=ny, spacing=spacing, margin=margin)
          md.Trajectory(sim.pos, sim.top, 0, sim.box, [90,90,90]).save(sim.pdb_cg)
          return

     # Each chain goes at the *centre* of its own cell (the +0.5), offset by the
     # margin — exactly what tools.new_simulation.lattice_positions previews.

     ibead = 0
     i = 0
     for comp in components:
          # CALVADOS builds each chain's starting conformation (build_compact)
          # as a serpentine through a cube of side ~cbrt(N) bonds, centred on
          # the origin — so bead 0 is the cube's lowest corner, at
          # -0.5*cbrt(N)*0.38 nm in x, y AND z (-1.9 nm for 1040 residues).
          # Until 2026-09 only the x/y centre of that cube was moved onto the
          # lattice point and z was shifted by a flat 2.0: the anchor ended up
          # ~2 nm off its lattice point sideways and at z = 0.1-1.4 nm, i.e.
          # BELOW the wall onset, permanently (it is fixed), with ~45 % of all
          # beads starting under the wall. Shifting by bead 0's own position
          # puts the anchor exactly on the lattice point at z_anchor and, since
          # bead 0 is the minimum corner, every other bead above it.
          xinit_anchored = comp.xinit - comp.xinit[0]

          for idx in range(comp.nmol):
               j = ibead + comp.nbeads

               x0 = margin + (i % nx + 0.5) * spacing
               y0 = margin + (i // nx + 0.5) * spacing

               # Bead 0 is the "Z"-tagged bead. CALVADOS gives "Z" a molecular
               # weight of -2 which the +2 N-terminus patch cancels to exactly
               # 0, and OpenMM holds zero-mass particles completely fixed — so
               # this lattice point is where the chain stays grafted for the
               # whole run, in x, y and z.
               pos = xinit_anchored + np.array([x0, y0, z_anchor])
               sim.pos[ibead:j] = pos

               ibead = j
               i += 1

     # The wall is stiff (wall_k), so a bead that starts under it would be
     # launched: refuse to write a starting structure with one. Cheap, and it
     # pins down the build_compact assumption above if CALVADOS ever changes it.
     lowest = float(np.asarray(sim.pos)[:, 2].min())
     if lowest < z_wall:
          raise RuntimeError(f"starting structure has a bead at z = {lowest:.3f} nm, below the "
                             f"wall onset at {z_wall} nm")

     md.Trajectory(sim.pos, sim.top, 0, sim.box, [90,90,90]).save(sim.pdb_cg)


# Job settings for Delft Blue
partition = "gpu-a100"
runtime = "1:30:00"
cpu_per_task = "18"

sim_name = Path(__file__).parent.name

box = [20.0, 20.0, 50.0]
N_save = 7000
N_frames = 5714

# OpenMM runs in nm. The surface is a one-sided harmonic wall on every bead,
# step(z_wall - z) * 0.5 * wall_k * (z_wall - z)^2: nothing below z_wall, and
# the "Z"-tagged residue 0 of every chain is held fixed at z_anchor, just above
# it. wall_k is what makes it impenetrable — see WALL_K in
# tools/new_simulation.py for the numbers behind it. metadata.csv records the
# full expression (ext_force_expr).
z_wall = 1.9
wall_k = 5000.0
z_anchor = 2.0

sequences: dict[str, str] = {
     "13r5g2d4": "LQLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGKGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPGIGVPGKGVPGMGVPGIGVPGIGVPGIGVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGKGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPGIGVPLDASTVYAVTGRGDSPASSAASAVPGIGVPGIGVPGIGVPGKGVPGIGVPGMGVPGIGVPGIGVPGIGVPGIGVPGIGVP"
}
if mode == "brush":
     # Residue 0 becomes the "Z" anchor bead (zero mass, held fixed by OpenMM).
     # In the lysine-attached modes it stays the residue the sequence says.
     sequences["13r5g2d4"] = f"Z{sequences['13r5g2d4'][1:]}"


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
          ext_force_expr = f'step({z_wall}-z)*0.5*{wall_k}*({z_wall}-z)^2',
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

     # Same for surface.yaml: its absence is what makes a run brush mode.
     write_surface_settings(runtime_dir,
                            SurfaceSettings.from_dict(surface) if surface else None)

     # Settings snapshot next to the trajectory-to-be. run.py rewrites it when
     # the run finishes (status/frame count); `sim metadata` refreshes it any time.
     write_metadata(runtime_dir)
