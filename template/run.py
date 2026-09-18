from calvados.sim import Sim
from yaml import safe_load
from ..prepare import build_sim # type: ignore
from tools.metadata import dcd_n_frames, write_metadata, write_timing
from datetime import datetime, timezone
from pathlib import Path
import time
import openmm

runtime_dir = Path(__file__).parent.resolve()

with open(runtime_dir / "config.yaml", 'r') as stream:
     config = safe_load(stream)

available_platforms = {openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())}
if config.get('platform') not in available_platforms:
     print(f"Platform '{config.get('platform')}' not available on this machine ({', '.join(sorted(available_platforms))}); falling back to CPU.")
     config['platform'] = 'CPU'

with open(runtime_dir / "components.yaml", 'r') as stream:
     components = safe_load(stream)

sim = Sim(runtime_dir, config, components)
sim.build_system()

build_sim(sim)

# How long the MD itself takes, wall-clock. Nothing CALVADOS writes records it,
# and it is what decides whether a batch of sequences is better sent as one
# long job or several short ones — so it is measured here and lands in
# metadata.csv (run_wall_seconds, ns_per_hour) with everything else about the
# run. Only simulate() is inside the clock: system building and, on DelftBlue,
# waiting in the queue are not part of what an extra sequence would cost.
started = datetime.now(timezone.utc)
clock = time.monotonic()

# Reactive crosslinking is opt-in. prepare.py writes runtime/crosslink.yaml only
# when it was asked for (tools.crosslink.CROSSLINK_FILENAME — spelled out here so
# a run with no crosslinking never even imports the module), and without that
# file this is the stock CALVADOS path, byte for byte as it was before the
# feature existed: no extra force, no random numbers drawn, same trajectory.
#
# Attaching chains through their lysines (tools.surface, modes "free" and
# "preattached") is the same kind of opt-in: prepare.py writes surface.yaml only
# for those modes, and without it the run is brush mode on the path above.
if (runtime_dir / "crosslink.yaml").is_file() or (runtime_dir / "surface.yaml").is_file():
     from tools.crosslink import load_settings, run_reactive
     from tools.surface import load_settings as load_surface_settings

     crosslink = load_settings(runtime_dir)
     surface = load_surface_settings(runtime_dir)
     if crosslink is None and surface is None:
          sim.simulate()
     else:
          run_reactive(sim, crosslink, runtime_dir, surface=surface)
else:
     sim.simulate()

elapsed = time.monotonic() - clock
write_timing(
     runtime_dir,
     seconds=elapsed,
     started=started.isoformat(timespec="seconds"),
     finished=datetime.now(timezone.utc).isoformat(timespec="seconds"),
     frames=dcd_n_frames(runtime_dir / f"{config['sysname']}.dcd"),
)

# The run is over, so metadata.csv can now record what actually came out of it
# (status, frames written, how long it took) on top of the settings prepare.py
# already stored.
write_metadata(runtime_dir)
