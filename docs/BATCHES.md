# Batches of runs from one CSV

`analyze.ipynb`'s first section plans a whole **set** of simulations at once:
one CSV row per run, one job per row. Batch CSVs live in `batches/`
(`batches/Example-csv.csv` is the template), and each holds only what differs
between runs:

```
Name ;Sequence;concentration;nmolecules;steps;crosslink distance;crosslink prob;mode;preattached fraction;walltime
example-v20-dilute;V20;0.02;16;7070000;;;;;0:40:00
example-block-xl;(GAGAGAGA) K1 V3 (VPGFG) V4;0.02;25;7070000;0.8;1.0;;;2:00:00
example-k4-free;K1 V4 K1 V4 K1 V4 K1 V4;0.02;16;7070000;;;free;;1:00:00
example-k4-pre;K1 V4 K1 V4 K1 V4 K1 V4;0.02;16;7070000;;;preattached;0.4;1:00:00
```

* **Any column may be missing, and any cell may be blank** — it then falls back
  to the value set once in the planning cell (`shared = Shared(...)`), which is
  also where everything that never varies within a batch lives: the margin, the
  length measure, the crosslink bond parameters, the SLURM partition.
* **`concentration` is always chains/nm²** here. The other two spacing modes
  (fraction of the chain's own length, grafted mass per area) are still in `sim
  new` and `create_simulation`, but a batch mixing modes wouldn't be comparable
  row to row.
* **`Sequence`** is either a raw single-letter sequence or the `tools/elibpy.py`
  shorthand: `V20` is 20 × `VPGVG`, `(VPGFG)` is inserted literally, so
  `K1 V3 (VPGFG) V4` is one lysine block, 3 × `VPGVG`, one `VPGFG`, 4 × `VPGVG`.
* **`walltime`** is per run (SLURM `HH:MM:SS`), so a long sequence and a short
  one can share a CSV without either waiting for the other's limit.
* **`mode`** is how the chains attach to the surface — blank or `brush` (the
  default, unchanged pipeline), `free` or `preattached` — see
  [Attaching through lysines](MODELS.md#attaching-through-lysines-free-and-preattached-modes).
  `preattached fraction` is only read for `preattached` rows (blank falls back to
  `surface_preattached_fraction` in `Shared`), so one CSV can sweep it.

The planning cell prints a comparison table of every run — box, spacing,
density, steps, frames, walltime — and, as 0/1 opt-ins, the full stats block
(`SHOW_STATS`) and the to-scale space-filling plot of the grafting lattice
(`SHOW_PLOTS`) for each one. Put a `0` behind either to skip it for a long CSV.

The table's footer adds up the walltimes two ways: the sum (what one job running
the rows back-to-back would have to ask for) and the maximum (the longest of
several separate jobs). Compare those against what earlier runs actually cost —
see `run_wall_seconds` / `ns_per_hour` under [Run metadata](SIMULATIONS.md#run-metadata).

## Names, and the list of used ones

`used-runs.csv` in the project root is the running list of every name and
sequence that has been submitted. The submit cell checks the whole CSV against
it **before preparing or queueing anything** — a name is refused if it repeats
inside the CSV, already exists as a folder under `simulations/`, or is already
in `used-runs.csv` — and appends the runs to it once `sbatch` has accepted them.
Nothing is recorded for a batch you only prepared, or for the local smoke-test
cell: the list is jobs that actually went out.

A reused *sequence* under a new name is reported as a note, not an error —
running the same ELP again at another concentration is a legitimate thing to do.

```bash
sim used              # show the list
sim used --backfill   # add any simulation on disk that isn't in it yet
```

`sim submit <name>` appends to the same list, so a one-off submission outside the
notebook is just as taken as one from a CSV.
