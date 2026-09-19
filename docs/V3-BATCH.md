# `block-runs-v3.csv` — the full 4 x 3 x 2 grid

24 runs: four block sequences, each in all three attachment modes, each with
lysine-lysine crosslinking on and off.

|  | brush | free | preattached |
| --- | --- | --- | --- |
| **K–K on** | `<family>-brush-xl-v3` | `<family>-free-xl-v3` | `<family>-pre-xl-v3` |
| **K–K off** | `<family>-brush-noxl-v3` | `<family>-free-noxl-v3` | `<family>-pre-noxl-v3` |

Families: `monoblock` (720 residues), `triblock-32` (880), `triblock-64` (1040),
`tetrablock` (1040). Sequences are copied verbatim from `block-runs.csv`, so
they are the same molecules as the v1 and v2 batches.

**`brush-noxl` is new.** The original `block-runs.csv` had only 20 rows because
that combination was never included, which left the brush column without its
control: there was no way to separate what crosslinking does from what the
brush geometry does. v3 closes the grid.

## What changed from v2

| setting | v1 / v2 | v3 | why |
| --- | --- | --- | --- |
| concentration | 0.02 chains/nm² | **0.015** | asked for |
| crosslink prob | 1.0 | **0.1** | asked for — see below |
| steps | 7,070,000 (70.7 ns) | **19,999,000 (200 ns)** | asked for |
| walltime | 3:00:00 | **4:00:00** | 2.83x the steps |
| crosslink distance | 0.8 nm | 0.8 nm | unchanged |
| nmol | 16 | 16 | unchanged |
| preattached fraction | 0.4 | 0.4 | unchanged |

Geometry that follows from the concentration: spacing 8.165 nm (was 7.071), box
32.66 x 32.66 nm (was 28.28). Box heights are unchanged at 85 / 105 / 120 nm,
since those scale with contour length, not density.

`steps` lands on 19,999,000 rather than 20,000,000 because `plan_steps` rounds
to a whole number of 7000-step frames: 2857 frames instead of v2's 1010.

## Walltime

v2 measured 13-18 min for 7,070,000 steps. Scaled by 2.83 that is **37-51 min**,
so 4 h leaves about 4.7x margin. Two of the changes push the real cost down
further: a lower grafting density means fewer neighbour pairs per bead, and
`prob = 0.1` means roughly ten times fewer bonds forming, hence far less of the
ramp sub-chunking that is the reactive loop's main overhead.

A job that somehow hits the limit resumes from its checkpoint with its network
intact, so a timeout costs a resubmission, not the run.

## Why `crosslink_prob = 0.1` is the most interesting change

With `prob = 1.0` every pair inside the capture radius reacted at the first
check, which is the diffusion-limited extreme. Real ELP gelation takes minutes
against nanoseconds of chain motion, i.e. it is strongly *reaction*-limited
(`docs/LITERATURE.md` LIT-5). Dropping to 0.1 moves the model toward that
regime.

It matters for more than realism. At `prob = 1.0` the reaction rate was set by
`crosslink_check_every`, a numerical parameter, so the result depended on how
often the code happened to look. With `prob < 1` the rate is carried by a
physical parameter instead, and the check interval stops being a hidden knob.

Expect **fewer crosslinks and a different intra/inter split**. The literature
(Zhong et al., Science 2016, and the PNAS follow-up) finds the primary-loop
fraction falling from about 35% to 19% on moving from fast to slow addition,
with a roughly sixfold effect on modulus. The v2 runs at `prob = 1.0` sat at
40-64% intra, so a drop is the expected direction and is the thing to look for.

## One thing to decide before submitting

`crosslink_start_step` and `surface_start_step` live in the notebook's `Shared`
block, not in the CSV, and are currently **1,000,000** — carried over from the
70 ns runs, where it was 14% of the trajectory. In a 200 ns run it is 5%.

The v2 equilibration analysis measured burn-in of up to **52 ns** (5.2M steps)
for these systems. So bonds formed from step 1,000,000 are still being made in
a brush that is far from relaxed. Since the whole point of 200 ns is better
statistics, consider raising it to around **5,000,000** so the network forms in
an equilibrated brush.

Left at 1,000,000 here because it was not part of the request, and because
changing it also changes what the run means.

## Still true, and still not fixed by any of this

At `crosslink_valence = 1` every lysine is a degree-≤2 node in the network
graph, so the crosslinks form paths and rings but never branch points. The
topology statistics are meaningful; the modulus is zero by construction until
f-functional junctions exist (`docs/GAPS.md` LIT-1,
`docs/DESIGN-crosslinker-models.md`).

And `G1 A1` in the sequences still expands to `VPGGGVPGAG…`, not a literal
`GAGA` block (`docs/BUGS.md` BUG-BATCH-2). Unchanged from v1 and v2, so the
three batches remain comparable.

## Submitting

All 24 names are checked free against `used-runs.csv` and `simulations/`. Point
the notebook's `RUNS_CSV` at `block-runs-v3.csv`, or prepare and submit from the
CLI.
