# What is unfinished on this branch

The audit ran six parallel auditors and was cut short by an API credit limit
partway through. Five of the six were killed mid-task. Their **code changes are
on disk and the test suite passes**, but several of them died before writing up
what they had done, and two died before finishing their job.

So: **a quiet area in `docs/BUGS.md` does not mean a clean area.** This file
says where to resume.

## Finished and trustworthy

* `docs/LITERATURE.md` (651 lines, real DOIs, with an explicit list of claims the
  author could *not* verify) and the `LIT-1…LIT-11` gap entries. This is the only
  agent that ran to completion.
* The crosslink/surface reactor audit: `tools/crosslink.py`, `tools/surface.py`,
  `tests/test_crosslink.py`, `tests/test_surface.py`. 47 tests pass in ~7 s with
  no GPU. Invariants covered are listed at the end of `docs/BUGS.md`.
* The wall and anchor fixes (`BUG-WALL-1`, `BUG-WALL-2`), with a measured
  before/after on `zz-audit-wall-before` / `zz-audit-wall-after`.
* `docs/SMOKE.md` inputs: eight end-to-end CPU runs under `simulations/zz-audit-*`
  with their logs in `docs/smoke-logs/`, including a two-leg restart.

## Unfinished, in priority order

1. **`tools/network.py` has no CLI command and no `docs/NETWORK.md`.**
   The module is complete (765 lines) and `tests/test_network.py` passes (14
   tests), but the agent was killed at the exact moment it went to add
   `sim network <sim>` to `tools/cli.py`. The module is importable and usable
   from Python today. To finish: append a `network` command to `tools/cli.py`
   just above the entry-point section, and write `docs/NETWORK.md` describing
   each metric and its citation. **This is the highest-value remaining work** —
   it is the only thing that turns crosslink events into a stiffness estimate,
   which is the project's actual deliverable.
2. **`tests/test_geometry.py` and `tests/test_box.py` were never written.**
   The wall/placement agent died before writing them. So the free/preattached
   hairpin construction (`build_chain`, `_hairpin_plan`, `check_geometry`) and
   the box planning (`plan_box`, `plan_lattice`, `plan_steps`,
   `mass_per_area`) have **no unit tests at all**. The intended cases are in
   that agent's brief: bond lengths exactly 0.38 nm, no non-bonded pair closer
   than that, pins at local z = 0, for n in {5, 10, 37, 100, 300, 720, 1040} and
   awkward pin sets (first, last, adjacent, all).
3. **`docs/SMOKE.md` was never written**, though every run it summarises exists.
   `docs/smoke-logs/analysis.log` and `analysis.json` hold the per-run numbers;
   the table just needs assembling. The per-invariant PASS/FAIL sweep reported
   `FAILS: none` for all eight runs.
4. **The notebook/batch/metadata audit stopped early.** It had read all 24
   notebook cells and fixed what is in `BUG-AN-*`, `BUG-BATCH-*`, `BUG-META-1`,
   but had not reached: the equilibration/burn-in statistics (N_eff,
   autocorrelation), whether Rg is computed on wrapped coordinates (it may be —
   this would be a real bug, and it is the observable the equilibration verdict
   rests on), RMSD alignment for grafted chains, and the DCD header parsing in
   `tools/metadata.py`. **Check Rg-on-wrapped-coordinates first.**
5. **Nothing was done about the crosslinker chemistry.** `docs/LITERATURE.md` §7
   and `docs/DESIGN-crosslinker-models.md` specify two models (implicit
   f-functional junctions; explicit crosslinker beads) and both agree the
   implicit one should be built first. Neither is implemented. This is the
   `LIT-1` blocker: with pairwise K–K bonds at any valence, the crosslinks are
   degree-≤2 graph nodes and cannot form a load-bearing network at all, so the
   modulus a finished `tools/network.py` computes would be zero by construction.

## The one thing to decide before resuming

`BUG-BATCH-2`: the 20 queued jobs (265624–265643) expand `G1 A1 G1 A1 …` to
`VPGGGVPGAG…`, not to the `GAGAGAGA` block the run names imply. That is a
scientific question about intent, not a code bug, and it is the user's call.
Everything else on this branch can proceed either way.

## How to resume safely

Read `docs/AGENT-RULES.md`. In short: work only in this worktree
(`/home/tnartey/ELP-Simulations-robust`, branch `robust-audit`), use
`/home/tnartey/ELP-Simulations/.venv/bin/python` with the worktree as cwd, never
`uv run`, never `sbatch`, and do not touch `/home/tnartey/ELP-Simulations` while
the batch is queued — those jobs import `tools/` from it at start-up.
