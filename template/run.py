from calvados.sim import Sim
from yaml import safe_load
from ..prepare import build_sim # type: ignore
from tools.metadata import write_metadata
from pathlib import Path
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

# Reactive crosslinking is opt-in. prepare.py writes runtime/crosslink.yaml only
# when it was asked for (tools.crosslink.CROSSLINK_FILENAME — spelled out here so
# a run with no crosslinking never even imports the module), and without that
# file this is the stock CALVADOS path, byte for byte as it was before the
# feature existed: no extra force, no random numbers drawn, same trajectory.
if (runtime_dir / "crosslink.yaml").is_file():
     from tools.crosslink import load_settings, run_reactive

     crosslink = load_settings(runtime_dir)
     if crosslink is None:
          sim.simulate()
     else:
          run_reactive(sim, crosslink, runtime_dir)
else:
     sim.simulate()

# The run is over, so metadata.csv can now record what actually came out of it
# (status, frames written) on top of the settings prepare.py already stored.
write_metadata(runtime_dir)
