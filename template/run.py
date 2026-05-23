from calvados.sim import Sim
from yaml import safe_load
from ..prepare import build_sim # type: ignore
from pathlib import Path

runtime_dir = Path(__file__).parent.resolve()

with open(runtime_dir / "config.yaml", 'r') as stream:
     config = safe_load(stream)

with open(runtime_dir / "components.yaml", 'r') as stream:
     components = safe_load(stream)

sim = Sim(runtime_dir, config, components)
sim.build_system()

build_sim(sim)

sim.simulate()