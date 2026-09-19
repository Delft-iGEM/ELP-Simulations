# GAPS

Branch robust-audit. Entries are appended by several auditors; ids are prefixed by area.


## LIT — modelling gaps from the literature review (see docs/LITERATURE.md)

* **LIT-1 (blocker). Pairwise K–K bonds cannot form a network at any `crosslink_valence`.**
  A lysine with valence 1 or 2 is a degree-≤2 graph node; degree-2 nodes are topologically
  invisible, so the crosslinks build paths and rings, never branch points, and contribute
  nothing to the modulus. THPP is f=3 (10.1021/bm061059m) and THPC is f=4
  (10.1021/bm3015279). Needs an f-functional junction (implicit is enough — see LITERATURE §7.2).
* **LIT-2. The surface wall is effectively transparent.** `step(z_wall-z)*0.5*(z_wall-z)^2` in
  an OpenMM CustomExternalForce is k = 1 kJ/mol/nm^2: 0.1 nm of penetration costs 0.002 kT and
  reaching z = 0 costs 0.74 kT. The box is periodic in z, so an escaping bead reappears at the
  top. k ≈ 2500 kJ/mol/nm^2 makes 0.1 nm cost 5 kT. Possibly a BUGS.md entry too.
* **LIT-3. Crosslink bond k = 2000 kJ/mol/nm^2 is ~100× too stiff** for what it represents
  (two lysine side chains + spacer, contour ≈1.5 nm). FJC entropic stiffness ≈20; k = 2000
  gives sigma = 0.035 nm, only 4× softer than CALVADOS's covalent backbone bond (8033).
  Recommend 250. r0 = 0.6 nm is, by contrast, correct (FJC predicts 0.61 nm).
* **LIT-4. Capture radius 0.8 nm is ~half the geometrically reachable CA–CA distance.**
  Lysine CA→NZ is 0.63 nm each side, so a THPP bridge spans ≈1.6 nm; Merkley et al. 2014
  (10.1002/pro.2458) show only ~19% of real crosslinks satisfy the naive static criterion and
  recommend hard-max × 1.1–1.25. Suggest 1.2 nm default, sweep 0.8/1.2/1.6.
* **LIT-5. React-on-first-contact is the wrong Damköhler limit.** Experimental ELP gelation is
  1–27 min (10.1021/bm3015279, 10.1021/bm061059m) against ns–µs chain relaxation, i.e. strongly
  reaction-limited; the code implements the diffusion-limited extreme. The kinetic regime
  changes primary-loop fraction 35%→19% and modulus ~6× (10.1073/pnas.1620985114).
  `tools/surface.py` already has a per-check `prob`; `tools/crosslink.py` does not — add it.
* **LIT-6. No modulus is ever computed.** The project's deliverable (stiffness vs sequence) is a
  pure post-processing job on data already written: build the crosslink graph, prune dangling
  ends, count nu/mu/xi, G_phantom = xi kT/V with V = footprint × brush height, minus the
  primary-loop correction (c_f ≈ 2.56 for f=3, 3.06 for f=4). Target range 0.25–46 kPa.
* **LIT-7. CALVADOS2 has never been validated on ELPs or on grafted/surface systems.** No
  CALVADOS–ELP paper found (searched PubMed + web, 2023–2026); recent extensions are PEG, PTMs
  and RNA. Every run is an extrapolation outside the validation set and should say so.
* **LIT-8. Langevin friction 0.01 ps^-1 is ~10^6× below water's Stokes friction** for a
  residue-sized bead — the dynamics are near-inertial. Event *times* are meaningless (already
  documented); less obviously, event *ordering* is only safe where it is set by distance and
  topology, not where it is set by how long a chain takes to travel.
* **LIT-9. Crosslinker stoichiometry is not a variable.** Experiment is non-monotonic in
  crosslinker:lysine (THPC: 0.5:1 → 250 Pa, 1:1 → 2200 Pa, >1:1 → slower gelation,
  10.1021/bm3015279) because excess crosslinker caps amines singly. No implicit model can
  reproduce this; it is the one thing an explicit-crosslinker-bead mode would buy.
* **LIT-10. pH trends will come out with the wrong sign.** Lim et al. measured |G*| rising
  5.8→9.8 kPa and 25.8→45.8 kPa from pH 7.5 to 12 (deprotonated Lys is the nucleophile). In
  CALVADOS, raising pH only discharges K and loosens the brush. Do not compare pH trends.
* **LIT-11. No LCST.** Already in the docstrings, and correct; the citable alternatives are
  HPS-T (10.1021/acscentsci.9b00102), the ELP-specific Baul/Dzubiella model
  (10.1021/acs.biomac.0c00546, whose authors aim it at "ELP networks and hydrogels"), or TEA
  bolted onto CALVADOS (10.1021/jacsau.6c00523, validated on 16 ELP sequences).

## RUN — performance and silent-failure gaps found from the 2026-09-19 batch

* **RUN-1. The reactive loop is 6–9.5x slower than the plain path.** Measured on
  the block batch: 7988 it/s non-reactive vs 1350 (monoblock) and 841
  (tetrablock) it/s reactive; 7.43e3 vs 1.02e3 ns/day. Cause: the loop calls
  `context.getState(getPositions=True)` every `check_every` = 1000 steps, which
  forces a full GPU→CPU sync of all positions, and `updateParametersInContext`
  re-uploads all 32 640 bond parameters on every change. Six of the twelve
  reactive jobs hit their 2 h walltime as a result. Either raise `check_every`
  (1000 steps = 10 ps is far finer than the chemistry justifies — see LIT-5),
  or restrict the pair list by a setup-time cutoff, or budget ~3 h.
* **RUN-2. Nothing detects a diverged simulation.** OpenMM does not error on
  coordinates of 1e13 nm; the run continues to the end, `metadata.csv` records
  "completed", and only the final PDB write fails — a run killed by walltime
  before that point reports no error at all. `report_potential_energy` is off,
  so the `StateDataReporter` log holds no energy trace to inspect afterwards
  either. **Add a cheap finiteness/extent check** at each checkpoint (max |z| vs
  box, or the potential energy) that fails the run loudly, and turn the energy
  reporting on. This is what let BUG-XL-7 pass unnoticed.
* **RUN-3. Smoke tests never produced an inter-chain crosslink.** Every audit
  smoke run and every earlier reactive test formed only intra-chain bonds in
  small boxes, so no bond ever spanned a periodic boundary and BUG-XL-7 was
  invisible to the whole test suite. Any future smoke matrix must include a
  case with chains close enough across a box edge to bond through it.
