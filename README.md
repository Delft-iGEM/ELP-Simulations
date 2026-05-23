# iGEM 2026 TU Delft

CALVADOS / OpenMM simulations of elastin-like polypeptides (ELPs).

## Setup

This project is managed with [`uv`](https://docs.astral.sh/uv/). Install it
by following the [official instructions](https://docs.astral.sh/uv/getting-started/installation/),
then clone the repo and run `uv sync`.

### Personal machine

```bash
git clone git@github.com:Delft-iGEM/ELP-Simulations.git && cd igem-2026
uv sync
```

In VSCode, select `.venv/bin/python` as your interpreter

For a terminal, activate the env with:

```bash
source .venv/bin/activate
```

### DelftBlue

Add the following to your `~/.bashrc` to redirect `uv`'s caches off the small
`/home` quota and install `uv` on first login:

```bash
export UV_CACHE_DIR=/scratch/$USER/.uv-cache
export UV_PYTHON_INSTALL_DIR=/scratch/$USER/.uv-python
mkdir -p "$UV_CACHE_DIR" "$UV_PYTHON_INSTALL_DIR"

export PATH="$HOME/.local/bin:$PATH"
command -v uv >/dev/null || curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then, **from a login node** (compute nodes have no internet):

```bash
cd /scratch/$USER
git clone <repo-url> && cd igem-2026
uv sync --frozen
```

## Writing simulations

The simulations use [CALVADOS](https://github.com/KULL-Centre/CALVADOS) on
top of OpenMM. If you haven't worked with CALVADOS before, skim the examples
in their repo first.

To create a new simulation, copy an existing folder under `simulations/`
and rename it. The folder name should match `sim_name` and be descriptive —
e.g. `i70-surfaceattached` for 70 `VPGIG` repeats, surface-attached. You can
use `tools.elibpy.humanize_seq` for naming inspiration.

Inside the simulation folder, `__init__.py` defines the usual CALVADOS
`Config` and `Components` setup; edit parameters there. ELP sequences can
be built with:

```python
from tools.elibpy import build_sequence_with_features
```

You don't have to use `elibpy` — plain strings in the `sequences` dict
work fine too.

The `runtime/` folder is generated and **not committed**. You should **NOT**
need to edit anything inside it by hand: `prepare.py` regenerates it
idempotently, so running `prepare.py` on any machine produces the exact same result.

## Preparing and running simulations

Generate the `runtime/` folder for a simulation with:

```bash
python -m simulations.elp_surfaceattached.prepare
```

This writes `runtime/run.py` along with all input files needed to launch the
simulation. Run it locally with:

```bash
python -m simulations.elp_surfaceattached.runtime.run
```

or submit it on DelftBlue with a SLURM script that calls the same `run.py`.
