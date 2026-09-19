# Re-running the 12 failed block runs — `block-runs-failed-v2.csv`

## Why these 12

The batch of 2026-09-19 (SLURM 265624–265643) split perfectly along one line:

| | outcome |
| --- | --- |
| 8 runs **without** crosslinking (`*-noxl`) | all COMPLETED in 12–16 min, usable |
| 12 runs **with** crosslinking (`*-xl`) | 5 FAILED, 6 TIMEOUT, 1 "completed" but corrupt |

Every crosslinking run diverged to coordinates of ~1e13 nm within two reaction
checks of `crosslink_start_step`. The cause is BUG-XL-7 in `docs/BUGS.md`: the
crosslink `HarmonicBondForce` was built without
`setUsesPeriodicBoundaryConditions(True)`.

The reaction test uses the minimum image, so a pair 0.8 nm apart *across* a box
face is bonded. Without the flag OpenMM then evaluated that bond on the raw
separation — a whole box length. Measured on a 28.284 nm box with the real
settings (k = 2000 kJ/mol/nm², r0 = 0.6 nm):

| | bond energy | force |
| --- | --- | --- |
| without the flag | 722,750 kJ/mol = 296,543 kT | 53,768 kJ/mol/nm |
| with the flag | 40 kJ/mol = 16 kT | 400 kJ/mol/nm |

That is the whole bug and the whole fix. It is one line, now covered by
`tests/test_crosslink.py::test_crosslink_force_uses_periodic_boundaries`, which
was verified to fail when the line is removed.

The `*-noxl` runs never formed a lysine–lysine bond, so nothing ever evaluated
that force, which is exactly why they are fine.

## What this CSV is

`block-runs-failed-v2.csv` is `block-runs.csv` filtered to the 12 crosslinking
rows. **Only two fields differ from the originals**, verified cell by cell:

* `Name` gains a `-v2` suffix, because the old names are burned in
  `used-runs.csv` and as folders under `simulations/`;
* `walltime` is raised.

Sequences, concentration, nmol, steps, crosslink distance, crosslink prob, mode
and preattached fraction are byte-identical to the first batch. This is a clean
re-run of the same science on fixed code, not a new experiment.

| runs | residues | walltime was | walltime now |
| --- | --- | --- | --- |
| `monoblock-*-xl-v2` | 720 | 2:00:00 | 4:00:00 |
| `triblock-32-*-xl-v2` | 880 | 2:00:00 | 5:00:00 |
| `triblock-64-*-xl-v2`, `tetrablock-*-xl-v2` | 1040 | 2:00:00 | 6:00:00 |

### Why so much more walltime

The reactive loop ran 6–9.5x slower than the plain path in the first batch:
1350 it/s (monoblock) and 841 it/s (tetrablock) against 7988 it/s for the
non-reactive runs. Six jobs hit the 2 h limit because of it.

Those rates are a **lower bound on the cost of a correct run**. They were
measured on systems that had already exploded, and a diverged system is cheap
to integrate: every bead is astronomically far from every other, so the
neighbour list is empty and the nonbonded terms cost nothing. A healthy
reactive run does the full nonbonded work *and* the reaction checks, so expect
it to be slower than 841 it/s, not faster. The walltimes above carry that
margin. A run that still hits the limit resumes from its checkpoint with its
network intact (`crosslink_state.json`), so a timeout costs a resubmission, not
the run.

The slowness itself is not fixed. `docs/GAPS.md` RUN-1 has the cause: the loop
calls `context.getState(getPositions=True)` every `check_every` = 1000 steps,
forcing a full GPU→CPU position sync, and `updateParametersInContext`
re-uploads all 32,640 bond parameters on every change. Raising `check_every`
is the cheapest lever — 1000 steps is 10 ps, far finer than the chemistry
justifies (see `docs/LITERATURE.md` LIT-5).

## Before you submit — two open decisions

Neither is a bug, and neither was changed on your behalf.

1. **`G1 A1` is not a `GAGA` block** (BUG-BATCH-2). In the shorthand, `X<n>`
   means n copies of the pentapeptide `VPGXG`, so the "triblock" mid-blocks
   expand to `VPGGGVPGAG…`. If you meant a literal alanine-rich elastic block,
   the cell has to read `(GAGAGAGA)`, as `Example-csv.csv` does. Changing it
   changes the sequences, so it is your call — and the 8 good `*-noxl` runs
   used the same expansion, so changing it here would break comparability with
   them.
2. **The crosslink parameters are unchanged from the first batch**: distance
   0.8 nm, k = 2000 kJ/mol/nm², p = 1.0, valence 1.
   `docs/LITERATURE.md` recommends 1.2 nm, k ≈ 250, and p < 1, with citations.
   Those are real improvements but they are a *different* experiment. Re-run
   as-is first so there is a clean baseline, then sweep.

Also worth knowing before drawing conclusions from the result: with
`crosslink_valence` 1 or 2 the crosslinks cannot form a load-bearing network at
all, because every lysine is a degree-≤2 graph node and those never branch
(`docs/GAPS.md` LIT-1). The runs will produce valid crosslink topology
statistics. They will not produce a non-zero modulus until f-functional
junctions exist (`docs/DESIGN-crosslinker-models.md`).

## How to submit

Point the notebook's `RUNS_CSV` at `block-runs-failed-v2.csv` and run the
planning and submit cells, or use the CLI per run. Keep
`crosslink_start_step` and `surface_start_step` at 1,000,000 as before. The
names have already been checked against `used-runs.csv` and `simulations/`, and
all 12 are free.
