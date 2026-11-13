# iGEM 2026 TU Delft — ELP Simulations

CALVADOS / OpenMM simulations of elastin-like polypeptides (ELPs).

## Setup

This project is managed with [`uv`](https://docs.astral.sh/uv/). Install it
by following the [official instructions](https://docs.astral.sh/uv/getting-started/installation/),
then clone the repo and run `uv sync`.

### Personal machine

```bash
git clone git@github.com:Delft-iGEM/ELP-Simulations.git
cd ELP-Simulations
uv sync
```

In VSCode, select `.venv/bin/python` as your interpreter. From a terminal,
activate the env with:

```bash
source .venv/bin/activate
```

### DelftBlue

Install `uv` once on a login node by following the
[official instructions](https://docs.astral.sh/uv/getting-started/installation/).

Then add the following to your `~/.bashrc` to redirect `uv`'s caches off the
small `/home` quota:

```bash
export UV_CACHE_DIR=/scratch/$USER/.uv-cache
export UV_PYTHON_INSTALL_DIR=/scratch/$USER/.uv-python
mkdir -p "$UV_CACHE_DIR" "$UV_PYTHON_INSTALL_DIR"
export PATH="$HOME/.local/bin:$PATH"
```

Reload your shell (`source ~/.bashrc`), then clone the repo to `/scratch`
and sync. This must be done **from a login node** — compute nodes have no
internet access:

```bash
cd /scratch/$USER
git clone git@github.com:Delft-iGEM/ELP-Simulations.git
cd ELP-Simulations
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

Inside the simulation folder, `prepare.py` defines the usual CALVADOS
`Config` and `Components` setup; edit parameters there. ELP sequences can
be built with:

```python
from tools.elibpy import build_sequence_with_features
```

You don't have to use `elibpy` — plain strings in the `sequences` dict
work fine too.

For anything that needs to happen on the OpenMM `System` itself (custom
restraints, extra forces, etc.), use the `build_sim` hook in `prepare.py`.
It runs after `sim.build_system()` and before `sim.simulate()`, so the
system exists and is ready to be modified.

The `runtime/` folder is generated and not committed. You **should NOT**
need to edit anything inside it by hand: `prepare.py` regenerates it
idempotently, so running prepare on any machine produces the exact same result.

## Preparing and running simulations

After `uv sync`, the `sim` CLI is available inside the venv. Activate it once:

```bash
source .venv/bin/activate   # Linux / macOS
.venv\Scripts\activate      # Windows (PowerShell)
```

Then use `sim` for all simulation tasks:

```bash
sim list                          # list available simulations
sim prepare i70-surfaceattached   # generate runtime/ folder
sim run     i70-surfaceattached   # run locally
sim submit  i70-surfaceattached   # sbatch on DelftBlue
```

### Getting the radius of gyration for a sequence

For a quick answer without writing a `prepare.py`, open the local web page:

```bash
sim ui
```

Enter a sequence as blocks (`I18 K17 I17 K18 RGD`, where `K70` means 70 `VPGKG`
repeats), a raw sequence, or FASTA. Set the conditions, pick a run length, and run.
The page plots the cumulative Rg (the mean of all kept frames up to each time) above
the per-frame Rg, reports the final value with its standard error, and lets you
download both. The chain is a single free molecule, and the first 10% of frames are
excluded as equilibration. The default Fast preset simulates 14 ns. Runs execute one
at a time and are saved under `runs/`.

To compare runs, tick two or more finished ones and press **Race selected**. Their
cumulative Rg curves grow together on a shared time axis, beside renders of each
chain moving along its trajectory. To give a run a picture, open it and paste with
Ctrl+V. The picture rides ahead of that run's line, and the winner's fills the screen
at the finish. **Hide Rg & score** hides every run's Rg and score until you're ready to race.

Every run uses the same conditions (20 °C, 0.13 M, pH 7.5, Fast). Each run needs a
target Rg, and its score is shown with its Rg:
900 − 100 × (nm from target), plus a variation bonus for different listed blocks
(ELP I/K/E/F/A/C/V/Q/W × 5/10/20, RGD, KR12, SAAP, IKVAV, XBBXBX; 7+ types = 50),
plus 50 when the sequence has the adhesion domain for the tissue the target falls in
(bone < 4 nm: RGD, skin 4–5 nm: IKVAV, brain ≥ 5 nm: XBBXBX). The block list is in `tools/rg.py`.

### Shell autocomplete

Install completion once (detects your shell automatically — bash, zsh, fish, or PowerShell):

```bash
sim --install-completion
```

Restart your shell (or `source ~/.bashrc`), then **Tab** after `sim prepare ` or
`sim run ` to autocomplete simulation names.

To see the raw completion script without installing it:

```bash
sim --show-completion
```

### Running without the venv active

Prefix any command with `uv run`:

```bash
uv run sim prepare i70-surfaceattached
uv run sim run     i70-surfaceattached
```
