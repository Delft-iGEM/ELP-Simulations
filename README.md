# iGEM 2026 TU Delft — ELP Simulations

CALVADOS / OpenMM simulations of elastin-like polypeptides (ELPs).

## Quick start

The project is managed with [`uv`](https://docs.astral.sh/uv/getting-started/installation/).

```bash
git clone git@github.com:Delft-iGEM/ELP-Simulations.git
cd ELP-Simulations
uv sync
source .venv/bin/activate
```

On DelftBlue, set up `uv`'s caches on `/scratch` first. See [docs/SETUP.md](docs/SETUP.md).

## Usage

```bash
sim list                          # list simulations
sim prepare i70-surfaceattached   # generate runtime/
sim run     i70-surfaceattached   # run locally
sim submit  i70-surfaceattached   # sbatch on DelftBlue
```

To plan, submit and analyse a whole batch of runs from one CSV, use
`analyze.ipynb`.

## Repository layout

| Path | Contents |
| --- | --- |
| `simulations/<name>/` | One folder per run: `prepare.py`, plus generated `analysis/` and `runtime/` (not committed) |
| `batches/` | Batch CSVs, one row per run (`Example-csv.csv` is the template) |
| `tools/` | The `sim` CLI, model code and analysis scripts |
| `template/` | `run.py` and `job.sh` copied into each run |
| `tests/` | `uv run pytest tests -q` |
| `docs/` | Documentation and development notes |
| `residues_CALVADOS2.csv` | CALVADOS2 residue parameters |
| `used-runs.csv` | Every run name and sequence already submitted |

## Documentation

- [Setup](docs/SETUP.md): installation, DelftBlue, where trajectory data goes, shell autocomplete
- [Writing and running simulations](docs/SIMULATIONS.md): `prepare.py`, the `sim` CLI, chain spacing, box size, run metadata
- [Batches](docs/BATCHES.md): planning many runs from one CSV, and the list of used names
- [Opt-in models](docs/MODELS.md): reactive crosslinking, CALVADOS-TEA, ionic strength, lysine attachment modes
- [Analysis](docs/ANALYSIS.md): the notebook's Analyze section, and the scripts behind each result series

## Authors and acknowledgment

- **Titus Nartey** calibrated the simulations, and made and analysed the runs
  for functional motif accessibility.
- **Christiaan Brunnings** created the runs for hydrogel stiffness modulation
  and the lysine placement module.
- **Joep Dols** provided the technical knowledge to simulate
  temperature-dependent behaviour.
- **Nina Sostaric** guided us on making figures and setting up these
  simulations.
- **Stefan van Loonen** helped us brainstorm through problems and provided
  helpful tools.
- **Bora Ciner** built the simulation pipeline and cli for working with CALVADOS.

Built on [CALVADOS](https://github.com/KULL-Centre/CALVADOS) (KULL Centre) and
[OpenMM](https://openmm.org).

## License

This repository is licensed under the
[Creative Commons Attribution 4.0 International License](LICENSE) (CC BY 4.0).