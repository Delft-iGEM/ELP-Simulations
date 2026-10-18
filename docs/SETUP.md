# Setup

Installing the environment, DelftBlue specifics, and where simulation data is stored.

This project is managed with [`uv`](https://docs.astral.sh/uv/). Install it
by following the [official instructions](https://docs.astral.sh/uv/getting-started/installation/),
then clone the repo and run `uv sync`.

## Personal machine

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

## DelftBlue

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

## Where simulation data goes

Trajectories are big — a few hundred MB per simulation — and DelftBlue caps
`/home` at 30 GB while `/scratch` gives you 5 TB. So on DelftBlue the
`runtime/` folder is a **symlink**: the data itself is written to

```
/scratch/$USER/elp-data/oefeningen/<sim_name>/
```

`sim prepare` and `sim new` create the link for you, and everything else —
`run.py`, `job.sh`, `analyze.ipynb` — keeps using the familiar
`simulations/<sim_name>/runtime/` path, so nothing changes in day-to-day use.
`sim clean` deletes the data on `/scratch` along with the link.

Set `ELP_DATA_DIR` to put the data somewhere else. On a machine with no
`/scratch/$USER` (your laptop), `runtime/` stays a plain folder inside the repo
and nothing is symlinked.

Simulations created before this was set up keep their data in the repo until
you move it:

```bash
sim migrate <simulation>   # or --all
```

## Shell autocomplete

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

## Running without the venv active

Prefix any command with `uv run`:

```bash
uv run sim prepare i70-surfaceattached
uv run sim run     i70-surfaceattached
```
