# Example simulations

Three example runs, taken from the full set of simulations in
[Delft-iGEM/ELP-Simulations](https://github.com/Delft-iGEM/ELP-Simulations)
(`simulations/`). The full set and `figures/` are left out of this submission
to keep it small. Every one of them can be regenerated from its `prepare.py`.

| Folder | Trajectory |
| --- | --- |
| `monoblock-free-noxl` | 1010 frames, 16 chains x 720 residues, 70.70 ns |
| `triblock-64-pre-noxl` | 1010 frames, 16 chains x 1040 residues, 70.70 ns |
| `tetrablock-free-noxl` | 1010 frames, 16 chains x 1040 residues, 70.70 ns |

Each folder contains only `prepare.py`, the CALVADOS `Config`/`Components`
setup that defines the run. Its outputs (`runtime/` with the trajectory, and
`analysis/` from `analyze.ipynb`) are never committed. Run `prepare.py`
(see `docs/SIMULATIONS.md`) to regenerate them.

To use an example with the `sim` tools, copy its folder into `simulations/`.
