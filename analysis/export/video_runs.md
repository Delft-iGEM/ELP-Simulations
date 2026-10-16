# SubM4 runs for the videos

Written by `analysis/scripts/08_video_runs.py` (repository HEAD `6725a2d`). Nothing was rendered.

Frozen definitions used for the contact files: lysine CA beads, **12 Å**, minimum image in x, y and z, **free lysines only** (neither lysine surface-pinned), intra-chain pairs ≥ 4 residues apart. Production starts at frame **1428** (100 ns burn-in); the metrics average frames 1428 onwards.

Files per run in `video_runs/`: `<Name>_contacts.csv` (columns `frame, type, chain_A, resid_A, chain_B, resid_B, index_A, index_B`; one row per frame and contact, all frames including burn-in; frame is the 0-based DCD frame index; resid is the residue number in `top.pdb`, 1-based within the chain; index is the 0-based atom index = VMD `index`), `<Name>_pins.csv` (`chain, resid, k_index, index`; K1-K4 counted from the N-terminus), and copies `<Name>_top.pdb`, `<Name>_surface_state.json`. The DCDs are not copied.

Units: DCD and PDB in Å; z is the surface normal (wall at z = 19 Å). Chains are written whole.

## 1. `sub4-anchor-BASE-CASE-L72-uniform-C0p04-I0p20-rep1`

Base case (BASE anchor, replicate 1).

C = 0.0400 chains/nm², L = 72 repeats, spacing 20.8 % (15 repeats), I = 0.200 M. Inter-chain free contacts 0.412/frame, intra 0.070/frame, bridge-to-loop ratio 5.93; 19 of 64 lysines pinned.

| | |
|---|---|
| DCD | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/runtime/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1.dcd` |
| (resolves to) | `/scratch/cjbrunings/elp-data/oefeningen/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1.dcd` |
| top.pdb | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/runtime/top.pdb` |
| surface_state.json | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/runtime/surface_state.json` |
| settings | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/runtime/config.yaml` (CALVADOS/OpenMM, nm); surface: `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/runtime/surface.yaml`; planner: `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-base-case-l72-uniform-c0p04-i0p20-rep1/prepare.py` |
| production starts at frame | 1428 |
| total frames | 5714 |
| contacts file check | every frame matches the extracted counts; production means 0.4120 inter / 0.0695 intra = `runs_metrics.csv` |

## 2. `sub4-anchor-A02-L64-S7-C0p01-I0p19-rep1`

The lowest grafting density of all finished runs (C = 0.0101 chains/nm²). Same simulated density as sub4-anchor-A02-L64-S7-C0p01-I0p19-rep2, sub4-anchor-A02-L64-S7-C0p01-I0p19-rep3, sub4-extra5-L87-S18-C0p01-I0p29-03; ties go to the lower replicate number, then the name.

C = 0.0101 chains/nm², L = 64 repeats, spacing 10.9 % (7 repeats), I = 0.190 M. Inter-chain free contacts 0.047/frame, intra 0.232/frame, bridge-to-loop ratio 0.20; 19 of 64 lysines pinned.

| | |
|---|---|
| DCD | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/runtime/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1.dcd` |
| (resolves to) | `/scratch/cjbrunings/elp-data/oefeningen/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1.dcd` |
| top.pdb | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/runtime/top.pdb` |
| surface_state.json | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/runtime/surface_state.json` |
| settings | `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/runtime/config.yaml` (CALVADOS/OpenMM, nm); surface: `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/runtime/surface.yaml`; planner: `/home/cjbrunings/ELP-Simulations/simulations/sub4-anchor-a02-l64-s7-c0p01-i0p19-rep1/prepare.py` |
| production starts at frame | 1428 |
| total frames | 5714 |
| contacts file check | every frame matches the extracted counts; production means 0.0474 inter / 0.2322 intra = `runs_metrics.csv` |

## 3. `sub4-extra4-L76-S4-C0p06-I0p35-28`

The highest grafting density of all finished runs (C = 0.0589 chains/nm²). Note: flagged for drift on the primary response (|z| = 4.6; kept in the analysis, README §5). Note: influential in the frozen model (Cook's D = 0.05, studentised residual 1.3; kept, README §5).

C = 0.0589 chains/nm², L = 76 repeats, spacing 5.3 % (4 repeats), I = 0.349 M. Inter-chain free contacts 1.361/frame, intra 0.609/frame, bridge-to-loop ratio 2.23; 19 of 64 lysines pinned.

| | |
|---|---|
| DCD | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l76-s4-c0p06-i0p35-28/runtime/sub4-extra4-l76-s4-c0p06-i0p35-28.dcd` |
| (resolves to) | `/scratch/cjbrunings/elp-data/oefeningen/sub4-extra4-l76-s4-c0p06-i0p35-28/sub4-extra4-l76-s4-c0p06-i0p35-28.dcd` |
| top.pdb | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l76-s4-c0p06-i0p35-28/runtime/top.pdb` |
| surface_state.json | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l76-s4-c0p06-i0p35-28/runtime/surface_state.json` |
| settings | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l76-s4-c0p06-i0p35-28/runtime/config.yaml` (CALVADOS/OpenMM, nm); surface: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l76-s4-c0p06-i0p35-28/runtime/surface.yaml`; planner: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l76-s4-c0p06-i0p35-28/prepare.py` |
| production starts at frame | 1428 |
| total frames | 5714 |
| contacts file check | every frame matches the extracted counts; production means 1.3609 inter / 0.6092 intra = `runs_metrics.csv` |

## 4. `sub4-extra3-L149-S37-C0p02-I0p24-09`

Longest chain combined with widest lysine spacing (largest coded L + coded S): L = 149 repeats, spacing 24.8 % (37 repeats). No run has both maxima: the longest chain is sub4-extra4-L150-S28-C0p06-I0p10-14 (L = 150, spacing 18.7 %), the widest spacing is sub4-extra3-L149-S37-C0p02-I0p24-09 (24.8 %).

C = 0.0231 chains/nm², L = 149 repeats, spacing 24.8 % (37 repeats), I = 0.245 M. Inter-chain free contacts 0.137/frame, intra 0.014/frame, bridge-to-loop ratio 10.14; 19 of 64 lysines pinned.

| | |
|---|---|
| DCD | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra3-l149-s37-c0p02-i0p24-09/runtime/sub4-extra3-l149-s37-c0p02-i0p24-09.dcd` |
| (resolves to) | `/scratch/cjbrunings/elp-data/oefeningen/sub4-extra3-l149-s37-c0p02-i0p24-09/sub4-extra3-l149-s37-c0p02-i0p24-09.dcd` |
| top.pdb | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra3-l149-s37-c0p02-i0p24-09/runtime/top.pdb` |
| surface_state.json | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra3-l149-s37-c0p02-i0p24-09/runtime/surface_state.json` |
| settings | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra3-l149-s37-c0p02-i0p24-09/runtime/config.yaml` (CALVADOS/OpenMM, nm); surface: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra3-l149-s37-c0p02-i0p24-09/runtime/surface.yaml`; planner: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra3-l149-s37-c0p02-i0p24-09/prepare.py` |
| production starts at frame | 1428 |
| total frames | 5714 |
| contacts file check | every frame matches the extracted counts; production means 0.1372 inter / 0.0135 intra = `runs_metrics.csv` |

## 5h. `sub4-extra4-L138-S22-C0p01-I0p08-16`

High bridge-to-loop ratio (2.40) at similar density to sub4-extra2-L58-S4-C0p01-I0p43-05 (ratio 0.06); densities 0.0123 vs 0.0128. Chosen as the pair with the largest ratio difference among runs whose densities differ by ≤ 5%.

C = 0.0123 chains/nm², L = 138 repeats, spacing 15.9 % (22 repeats), I = 0.076 M. Inter-chain free contacts 0.080/frame, intra 0.033/frame, bridge-to-loop ratio 2.40; 19 of 64 lysines pinned.

| | |
|---|---|
| DCD | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l138-s22-c0p01-i0p08-16/runtime/sub4-extra4-l138-s22-c0p01-i0p08-16.dcd` |
| (resolves to) | `/scratch/cjbrunings/elp-data/oefeningen/sub4-extra4-l138-s22-c0p01-i0p08-16/sub4-extra4-l138-s22-c0p01-i0p08-16.dcd` |
| top.pdb | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l138-s22-c0p01-i0p08-16/runtime/top.pdb` |
| surface_state.json | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l138-s22-c0p01-i0p08-16/runtime/surface_state.json` |
| settings | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l138-s22-c0p01-i0p08-16/runtime/config.yaml` (CALVADOS/OpenMM, nm); surface: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l138-s22-c0p01-i0p08-16/runtime/surface.yaml`; planner: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra4-l138-s22-c0p01-i0p08-16/prepare.py` |
| production starts at frame | 1428 |
| total frames | 5714 |
| contacts file check | every frame matches the extracted counts; production means 0.0800 inter / 0.0334 intra = `runs_metrics.csv` |

## 5l. `sub4-extra2-L58-S4-C0p01-I0p43-05`

Low bridge-to-loop ratio (0.06) at similar density to sub4-extra4-L138-S22-C0p01-I0p08-16 (ratio 2.40); densities 0.0128 vs 0.0123. Chosen as the pair with the largest ratio difference among runs whose densities differ by ≤ 5%. Note: influential in the frozen model (Cook's D = 0.95, studentised residual -6.2; kept, README §5).

C = 0.0128 chains/nm², L = 58 repeats, spacing 6.9 % (4 repeats), I = 0.427 M. Inter-chain free contacts 0.040/frame, intra 0.640/frame, bridge-to-loop ratio 0.06; 19 of 64 lysines pinned.

| | |
|---|---|
| DCD | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra2-l58-s4-c0p01-i0p43-05/runtime/sub4-extra2-l58-s4-c0p01-i0p43-05.dcd` |
| (resolves to) | `/scratch/cjbrunings/elp-data/oefeningen/sub4-extra2-l58-s4-c0p01-i0p43-05/sub4-extra2-l58-s4-c0p01-i0p43-05.dcd` |
| top.pdb | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra2-l58-s4-c0p01-i0p43-05/runtime/top.pdb` |
| surface_state.json | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra2-l58-s4-c0p01-i0p43-05/runtime/surface_state.json` |
| settings | `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra2-l58-s4-c0p01-i0p43-05/runtime/config.yaml` (CALVADOS/OpenMM, nm); surface: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra2-l58-s4-c0p01-i0p43-05/runtime/surface.yaml`; planner: `/home/cjbrunings/ELP-Simulations/simulations/sub4-extra2-l58-s4-c0p01-i0p43-05/prepare.py` |
| production starts at frame | 1428 |
| total frames | 5714 |
| contacts file check | every frame matches the extracted counts; production means 0.0397 inter / 0.6400 intra = `runs_metrics.csv` |

## 6. Wet-lab constructs

**None selected.** No run of Monoblock, Triblock-32, Triblock-64 or Tetrablock is in `runs_metrics.csv`: every SubM4 run is a (VPGIG)ₙ chain with four VPGKG.

Folders with these construct names do exist from earlier submodules (140 under `simulations/`: monoblock 35, triblock-32 35, triblock-64 35, tetrablock 35), but 0 of them have a readable `runtime/` (trajectory): they hold only `prepare.py` and `analysis/` summaries. They were run under other protocols (box, length, burn-in) and are not part of this dataset, so nothing was chosen for them.
