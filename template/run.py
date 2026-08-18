from calvados.sim import Sim
from yaml import safe_load
from ..prepare import build_sim # type: ignore
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

sim.simulate()