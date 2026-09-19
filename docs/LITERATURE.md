# LITERATURE — parameter memo for the ELP crosslinking / surface pipeline

Branch `robust-audit`. Written to answer seven parameter questions with primary sources.
Scope rule (docs/AGENT-RULES.md): this file and appended entries in `docs/GAPS.md` are the
only things written; no code was touched.

**Conventions.** Distances in nm (the code's units), energies in kJ/mol.
kT at T = 293.15 K is **2.437 kJ/mol**. A harmonic bond `0.5 k (r-r0)^2` has
RMS fluctuation `sigma = sqrt(kT/k)`.

**Verification policy.** Every number below is either (a) quoted from a source with a DOI,
(b) an arithmetic consequence of such numbers, labelled *[derived]*, or (c) explicitly
flagged *[UNVERIFIED]* / *[inferred]*. Bibliographic metadata for the PubMed-indexed items
was retrieved from PubMed; DOIs are given as links. Nothing here is a citation I could not
resolve to a real record.

---

## 0. What the current model does (read from the code, for reference)

| thing | value in repo | where |
| --- | --- | --- |
| force field | CALVADOS2, 1 bead/residue, Ashbaugh–Hatch + Debye–Hückel | `.venv/.../calvados` |
| backbone bond | harmonic, r0 = 0.38 nm, k = 8033 kJ/mol/nm^2 | CALVADOS default |
| integrator | Langevin, dt = 10 fs, friction 0.01 ps^-1 | CALVADOS default |
| T / ionic / pH | 293.15 K / 0.19 M / 7.5 | `simulations/*/prepare.py` |
| wall | `step(z_wall-z)*0.5*(z_wall-z)^2`, z_wall = 1.9 nm | `prepare.py`, `Config(ext_force_expr=...)` |
| wall stiffness | **k = 1 kJ/mol/nm^2** (no prefactor in the expression) | same |
| crosslink capture radius | 0.8 nm CA–CA | `tools/crosslink.py:137-139` |
| crosslink bond | r0 = 0.6 nm, k = 2000 kJ/mol/nm^2, ramped over 500 steps | same |
| valence | 1 (bifunctional) or 2 per lysine | same |
| grafting density | 0.02 chains/nm^2 | `tools/new_simulation.py` |

---

## 1. Crosslinker geometry: functionality f, and the CA–CA distance a bond should represent

### 1.1 The reagents

| reagent | what it is | functionality f (amines bridged) | source |
| --- | --- | --- | --- |
| **THPP** | β-[tris(hydroxymethyl)phosphino]propionic acid — a **phosphine** with three hydroxymethyl arms; reacts with primary/secondary amines by Mannich-type condensation forming **P–CH2–N** links, water the only by-product | **3** | Lim, Nettles, Setton & Chilkoti 2007, [10.1021/bm061059m](https://doi.org/10.1021/bm061059m) |
| **THPC** | tetrakis(hydroxymethyl)phosphonium **chloride**; explicitly described as "tetra-functional", four hydroxymethyl arms; proceeds via a formaldehyde intermediate (2.3–9.2 mM free HCHO measured at encapsulation doses) | **4** | Chung, Lampe & Heilshorn 2012, [10.1021/bm3015279](https://doi.org/10.1021/bm3015279) |
| **glutaraldehyde** | OHC-(CH2)3-CHO, but in water exists in ≥13 interconverting forms (monomer, hydrate, oligomers); the reactive species is disputed and crosslinking is "largely developed through empirical observation" | nominally **2**, effectively **2–many** (oligomeric bridges) | Migneault *et al.* 2004, [10.2144/04375RV01](https://doi.org/10.2144/04375RV01) |
| **genipin** | iridoid; two distinct reactions — fast nucleophilic attack of an amine on the dihydropyran ring, slow amide formation at the ester — so each genipin joins **two** amines; further complicated by oxygen-radical **oligomerisation** of genipin (the blue colour), which makes longer multi-amine bridges | **2** nominal, longer/oligomeric in practice | Butler, Ng & Pudney 2003, [10.1002/pola.10960](https://doi.org/10.1002/pola.10960) |
| **HMDI** (hexamethylene diisocyanate) | used on elastin-/fibronectin-derived aECM protein films, crosslinked in DMSO, tensile moduli spanning native elastin | **2** | Nowatzki & Tirrell 2004, [10.1016/S0142-9612(03)00635-5](https://doi.org/10.1016/S0142-9612\(03\)00635-5) |
| **transglutaminase** | enzyme; forms the **zero-length** Nε-(γ-glutamyl)lysine isopeptide bond between a Gln side chain and a Lys side chain | **2**, and it is **K–Q**, not K–K | used for ELP gels by McHale, Setton & Chilkoti 2005, [10.1089/ten.2005.11.1768](https://doi.org/10.1089/ten.2005.11.1768) |
| **DSS / BS3** (reference only) | NHS esters, spacer arm **1.14 nm** fully extended | **2** | Merkley *et al.* 2014, [10.1002/pro.2458](https://doi.org/10.1002/pro.2458) |

> **The single most important line in this whole memo:** THPP is **f = 3** and THPC is **f = 4**.
> Neither is a bifunctional linker, and neither is what `crosslink_valence` currently models
> (see §2).

### 1.2 Translating a spacer length into a CA–CA distance

Merkley *et al.* 2014 ([10.1002/pro.2458](https://doi.org/10.1002/pro.2458)) is the reference
calculation. Their numbers, quoted:

* DSS/BS3 linker arm **11.4 Å** fully extended.
* "the summed lengths of the linker and the two lysine side chains **[24 Å]** represent a hard
  upper limit to the Cα–Cα distance in the crosslinked conformation" — i.e. **2 × ≈6.3 Å of
  lysine side chain** (CA→NZ), which is the 0.63 nm the task statement assumes. Good.
* But static structures badly under-predict what reacts: "only **∼19%** of observed crosslinks
  had a Nζ–Nζ distance shorter than the crosslinker arm of BS3 or DSS (11.4 Å)".
* From MD (the Dynameomics database, 807 proteins) they conclude a practical restraint of
  **26–30 Å between Cα atoms**, and note pairs starting "as great as ∼35–40 Å could still
  become crosslinked".

So the empirical rule is: **usable CA–CA distance ≈ (hard geometric maximum) × 1.1–1.25**, the
excess coming from side-chain rotamer flips and backbone motion, not from the linker.

Applying that rule to our reagents *[derived]*:

| linker | bridge length between the two Nζ | hard max CA–CA = 0.63 + 0.63 + bridge | Merkley-style working capture radius |
| --- | --- | --- | --- |
| transglutaminase (zero-length K–Q) | ~0.15 nm (one C–N bond) | ~1.4 nm | ~1.5 nm |
| **THPP / THPC** (N–CH2–P–CH2–N) | **~0.3 nm** *[derived from standard P–C 0.187 nm, C–N 0.147 nm bond lengths — UNVERIFIED against a crystal structure of the adduct]* | **~1.6 nm** | **~1.6–1.9 nm** |
| glutaraldehyde monomer (C5 bridge) | ~0.6 nm | ~1.9 nm | ~2.1–2.3 nm |
| genipin | ~0.7–0.9 nm *[UNVERIFIED — I did not find a measured Nζ–Nζ span for the genipin adduct]* | ~2.1 nm | ~2.4 nm |
| DSS/BS3 | 1.14 nm (measured) | 2.4 nm (their number) | **2.6–3.0 nm (their number)** |

**The current 0.8 nm capture radius is therefore about half the geometrically reachable
distance for a THPP crosslink.** It is not wrong as a *conservative* criterion — but it should
be understood as "the two CA beads are nearly touching", which is a much stricter condition
than "these two lysines can be bridged by THPP".

### 1.3 What bond stiffness represents a flexible spacer?

The crosslink in the CG model has to stand in for *two lysine side chains plus the linker* —
i.e. a short flexible tether of contour length Lc ≈ 1.5 nm (THPP) between the two CA beads.
Treat it as a freely jointed chain of Kuhn length b ≈ 0.25 nm (roughly two single bonds per
Kuhn segment). Then *[derived, standard FJC/Gaussian result]*:

```
<r^2> = Lc * b                -> r_rms = sqrt(1.5 * 0.25) = 0.61 nm
F(r)  = 3 kT r^2 / (2 Lc b)   -> k_entropic = 3 kT / (Lc b) = 3*2.437/0.375 = 19.5 kJ/mol/nm^2
```

Two conclusions:

1. **r0 = 0.6 nm is well chosen.** The FJC most-probable CA–CA separation for a THPP-length
   tether is 0.61 nm. Keep it. (For a longer linker — glutaraldehyde, genipin, BS3 — r0 should
   scale as sqrt(Lc·b): ≈0.75 nm, ≈0.8 nm, ≈0.9 nm respectively *[derived]*.)
2. **k = 2000 kJ/mol/nm^2 is ~100× too stiff** for what it represents. It gives
   sigma = sqrt(2.437/2000) = **0.035 nm** *[derived]* — a crosslink that barely breathes,
   only 4× softer than CALVADOS's own covalent backbone bond (8033 kJ/mol/nm^2,
   sigma = 0.017 nm; [10.48550/arXiv.2504.10408](https://doi.org/10.48550/arXiv.2504.10408)).
   Physically the crosslink is a *floppier* object than a peptide bond, not a comparable one.

A pure entropic spring (k ≈ 20) is also not right for a simulation bond, because a harmonic
form with that stiffness lets r wander to 1.5 nm and beyond for ~2 kT, where the real tether
is inextensible. The practical compromise is to pick k from a target fluctuation amplitude of
0.1–0.15 nm, which is the scale over which a real 1.5-nm tether is genuinely soft:

```
k = kT / sigma^2 :  sigma = 0.15 nm -> k = 108 ;  sigma = 0.10 nm -> k = 244 kJ/mol/nm^2
```

**Recommendation: k = 250 kJ/mol/nm^2** (sigma ≈ 0.10 nm), with 100 and 500 as the sweep
bracket. This is soft enough that a forming crosslink does not yank two chains together with
covalent force, and stiff enough to stay a topological constraint.

*Practical note for the ramp:* a softer k also reduces the energy injected when a bond switches
on, which is exactly what `crosslink_ramp_steps` exists to manage
(`tools/crosslink.py` docstring). Softer k makes the ramp less critical, not more.

---

## 2. Functionality and topology: why `valence = 2` is not `f = 3`

### 2.1 The topological distinction

`crosslink_valence` is the number of crosslink bonds **one lysine bead** may carry. Setting it
to 2 makes a lysine a **degree-2 node** — it can be a waypoint on a path, chaining
K–K–K in a line. A **linear chain of pairwise bonds carries no elasticity beyond the chains it
already connects**: a degree-2 node is topologically invisible, it can be smoothed away.

An f = 3 crosslinker is a **degree-3 node**: a single molecule (THPP) bonded to three
*different* lysines at once, forming a star junction. Only nodes of degree ≥ 3 create network
branch points, and only branch points contribute to the modulus.

Concretely, with N lysines:
* valence-1 pairwise bonds: builds a graph of degree ≤1 nodes → at best a linear multiblock.
  This models glutaraldehyde at 1:1, correctly.
* valence-2 pairwise bonds: a graph of degree ≤2 → still only paths and rings. **Still no
  branch points, still no gel.** Doubling the valence buys longer paths and *primary loops*,
  not a network.
* an explicit or implicit f = 3 junction: degree-3 nodes → percolation and a finite modulus.

So the current code, at either valence setting, **cannot form a load-bearing 3D network from
K–K bonds alone.** What holds the coating together in the runs is the surface anchoring plus
the chain connectivity, not the crosslinks. This is the biggest single modelling gap found
(logged in GAPS.md).

There is direct experimental evidence that this distinction is not cosmetic. Kawamoto, Zhong,
Wang, Olsen & Johnson (Macromolecules 2015, 48, 8980–8988,
[10.1021/acs.macromol.5b02243](https://doi.org/10.1021/acs.macromol.5b02243)) compared A2+B3
against A2+B4 click hydrogels; higher junction functionality is markedly more prone to cyclic
defects. Wang, Alexander-Katz, Johnson & Olsen (PRL 2016, 116, 188302,
[10.1103/PhysRevLett.116.188302](https://doi.org/10.1103/PhysRevLett.116.188302)) then showed
that *all* cyclic topologies collapse onto a **universal function of a single dimensionless
parameter** characterising the formation conditions, and that the whole topology is fixed by
the primary-loop fraction alone — "the one-to-one correspondence between the network topology
and primary loop fraction demonstrates that the entire network topology is characterized by
measurement of just primary loops".

### 2.2 Formulas a code implementation should use

**Step 1 — build the graph.** Nodes = grafting anchors, crosslink junctions and chain ends.
Edges = (a) chain strands between consecutive junctions on the same chain, (b) crosslink bonds,
(c) surface tethers. Store edge lengths in residues.

**Step 2 — prune.** Iteratively delete every node of degree 1 together with its edge
(this removes dangling ends and whole dangling branches). Keep only the component(s)
percolating across the periodic x,y box and/or connected to the surface. This is the
**elastically active network**; everything removed is sol + dangling material, which carries
no equilibrium stress.

**Step 3 — count.** On the pruned graph, for each connected component:

```
nu  = number of elastically active strands   (edges)
mu  = number of elastically active junctions (nodes)
xi  = nu - mu + 1          # cycle rank, per component  (= E - V + C overall)
```

**Step 4 — moduli.**

```
Affine  (upper bound):   G_aff = nu * kT / V
Phantom (lower bound):   G_ph  = (nu - mu) * kT / V = xi * kT / V
```
The phantom result is Flory's ([Proc. R. Soc. Lond. A 351, 351 (1976),
10.1098/rspa.1976.0146](https://doi.org/10.1098/rspa.1976.0146)); for a perfect f-functional
network mu = 2nu/f, so `G_ph = nu kT (1 - 2/f) / V` — **zero at f = 2**, which is the algebraic
statement of §2.1. For f = 3 the phantom modulus is 1/3 of affine, for f = 4 it is 1/2.

**Step 5 — loop correction.** Count primary loops L1 (a crosslink joining two lysines on the
same chain with no intervening junction, and any pair of parallel crosslinks between the same
two strands). Two equivalent-in-spirit corrections from the literature:

*Real Elastic Network Theory* (Zhong, Wang, Kawamoto, Olsen & Johnson, Science 2016, 353,
1264–1268, [10.1126/science.aag0184](https://doi.org/10.1126/science.aag0184)), as written in
the follow-up PNAS paper on loop-defect control
([10.1073/pnas.1620985114](https://doi.org/10.1073/pnas.1620985114)) for A2+B4:

```
G = G0 * (1 - 8/3 * phi_1)
```
with phi_1 the primary-loop fraction and G0 the defect-free modulus.

*Phantom-network form with an explicit loop term* (Panyukov-style derivation reviewed in
["On the Elasticity of Polymer Model Networks Containing Finite Loops",
arXiv:2103.16181](https://arxiv.org/abs/2103.16181)):

```
G_ph ≈ (xi - c_f * L1) * kT / V ,  c_f ≈ 2.56 (f=3), 3.06 (f=4), -> ≈4.2 as f -> inf
```

Use the second form: it plugs straight into the graph counts of Step 3 and is
functionality-aware.

**Step 6 — sanity/scale.** For a *surface coating* the relevant experimental observable is
usually the layer's effective shear modulus; V should be the volume actually occupied by the
brush (box footprint × measured brush height from `z_distribution.py`), not the whole
(mostly empty) simulation box. Getting this wrong changes G by the ratio of box height to
brush height, which in the current setups is a factor of 3–5.

**Miller–Macosko** (Macromolecules 1976, 9, 206–211,
[10.1021/ma60050a004](https://doi.org/10.1021/ma60050a004)) gives the mean-field version of the
same quantities, for cross-checking the graph counts against conversion p. For self-reacting
f-functional sites, the probability that looking out from a reacted group one finds a *finite*
chain, P, solves

```
P = (1 - p) + p * P^(f-1)
```
and a junction is elastically effective when at least 3 of its f arms lead to the infinite
network; the gel point is at `p_gel = 1/(f-1)` (= 0.5 for f = 3, 1/3 for f = 4). Mean-field
theory ignores loops by construction, which is precisely why the measured modulus falls below
it and why Step 5 exists.

---

## 3. Experimental calibration points: ELP gel moduli

Six datapoints, all from primary sources. Where a paper's own notation had to be interpreted
to get a lysine fraction, that is marked *[inferred]*.

| # | system / sequence | K content | [ELP] | crosslinker & ratio | modulus | source |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | ELP, VPGXG with X = Lys every **7** pentapeptides (else Val) | high | varied | **TSAT** (tris-succinimidyl aminotriacetate, f = 3) | **0.24–3.7 kPa at 7 °C**, **1.6–15 kPa at 37 °C** across the whole series (conc., Lys content, MW) | Trabbic-Carlson, Setton & Chilkoti 2003, [10.1021/bm025671z](https://doi.org/10.1021/bm025671z) |
| 2 | same study, X = Lys every **17** pentapeptides | low | varied | TSAT | lower end of the same ranges | ibid. |
| 3 | **ELP[KV7F-72]** | 1 K per 9 pentapeptides ≈ **2.2 mol%** *[inferred from the [K V7 F] repeat notation]* | ≈170 mg/mL final | **THPP**, pH 7.5 | **\|G*\| ≈ 5.8 kPa** (10 rad/s) | Lim, Nettles, Setton & Chilkoti 2007, [10.1021/bm061059m](https://doi.org/10.1021/bm061059m) |
| 4 | **ELP[KV2F-64]** | 1 K per 4 pentapeptides ≈ **5 mol%** *[inferred, same notation]* | ≈170 mg/mL | THPP, pH 7.5 | **\|G*\| ≈ 25.8 kPa**, i.e. **3.6–4.8× the KV7F gel**; rising to 45.8 kPa at pH 12 | ibid. |
| 5 | recombinant **elastin-like protein**, 5 wt% | designed Lys sites | **5 wt%** | **THPC** 1:1 amine ratio (gel time 6.7 ± 0.2 min; 27 ± 1.2 min off-ratio) | **G′ ≈ 2200 Pa**; at 0.5:1 only ≈250 Pa. Concentration series: 3 wt% **300 ± 22 Pa**, 5 wt% **1320 ± 23 Pa**, 10 wt% **3090 ± 16 Pa** | Chung, Lampe & Heilshorn 2012, [10.1021/bm3015279](https://doi.org/10.1021/bm3015279) |
| 6 | ELP block copolymers (monoblock ELP[KV7F-144]; ABA and BABA with non-crosslinkable ELP[VG7A8] middle blocks) | A block = 1 K per 9 pentapeptides | **200 mg/mL** | **THPP**, 1:1 HMP:Lys, gel time **< 1 min** | G′ rose **4000 → 6000 Pa** over 30 min of crosslinking; monoblock stiffest, long-B-block triblock softest | Lim, Nettles, Setton & Chilkoti 2008, [10.1021/bm7007982](https://doi.org/10.1021/bm7007982) |
| 7 | ELP + transglutaminase (K–Q, enzymatic) | — | — | tissue transglutaminase | dynamic shear modulus **0.28 → 1.7 kPa** over 4 weeks of chondrocyte culture (the rise is deposited cartilage matrix, not crosslinking) | McHale, Setton & Chilkoti 2005, [10.1089/ten.2005.11.1768](https://doi.org/10.1089/ten.2005.11.1768) |

**Genipin + ELP:** I could **not** verify a clean genipin-crosslinked-ELP modulus datapoint from
a primary source. Secondary sources point at "Thermosensitive elastin-derived polypeptide
hydrogels crosslinked by genipin", Int. J. Polym. Mater. 66(8) 2017,
[10.1080/00914037.2016.1217534](https://doi.org/10.1080/00914037.2016.1217534), with a
compressive modulus reportedly up to ≈39 kPa, but that is behind a paywall and I did not read
the paper — **treat as UNVERIFIED**. The genipin ELP datapoints circulating in reviews trace
back to acid-degraded elastin peptides rather than recombinant VPGXG ELPs, which is a different
material.

**What the table is good for.** Two robust trends to reproduce qualitatively:
1. **Modulus scales strongly with lysine density** — 2.2 → 5 mol% K gave a 3.6–4.8× stiffness
   jump at fixed concentration and crosslinker ratio (rows 3 vs 4). This is the single best
   calibration target for "different ELP sequences → different stiffness", which is the project
   goal.
2. **Modulus is non-monotonic in crosslinker:lysine ratio** — THPC at 1:1 gives ≈2200 Pa, at
   0.5:1 only ≈250 Pa, and going *above* 1:1 slows gelation again (row 5). Excess crosslinker
   caps amines singly instead of bridging them. A model with no explicit crosslinker species
   **cannot reproduce this at all** (see §7).

The absolute range to aim for is **0.25–46 kPa**, i.e. ~10^2.5–10^4.5 Pa.

---

## 4. Force-field validity

### 4.1 CALVADOS is the right family; it has no LCST

CALVADOS2 (Tesei & Lindorff-Larsen, Open Res. Europe 2:94,
[10.12688/openreseurope.14967.2](https://doi.org/10.12688/openreseurope.14967.2); original
parameterisation Tesei, Schulze, Crehuet & Lindorff-Larsen, PNAS 2021,
[10.1073/pnas.2111696118](https://doi.org/10.1073/pnas.2111696118)) is optimised against
single-chain dimensions of 55 IDPs plus 70 hydrophobicity scales, and validated for chain
compaction and phase-separation propensity "at different temperatures and salt concentrations".
CALVADOS 3 extends this to multi-domain proteins (Cao, von Bülow, Tesei & Lindorff-Larsen,
Protein Sci. 2024, [10.1002/pro.5172](https://doi.org/10.1002/pro.5172)).

**But the residue "stickiness" λ is temperature-independent.** Temperature enters only through
kT and through the Debye screening length. An ELP's defining behaviour — the LCST/inverse
transition, where hydrophobic contacts get *stronger* on heating — is exactly the part that a
fixed-λ model cannot produce. The `tools/crosslink.py` docstring already says this; it is
correct, and the citations for it are:

* **Dignon, Zheng, Kim & Mittal, ACS Cent. Sci. 2019, 5, 821–830,
  [10.1021/acscentsci.9b00102](https://doi.org/10.1021/acscentsci.9b00102)** — the HPS-T model:
  a knowledge-based amino-acid potential with an explicitly **temperature-dependent**
  solvent-mediated interaction, which distinguishes >35 UCST from LCST IDPs. This is the
  canonical "if the transition matters, use this" reference.
* **Baul, Bley & Dzubiella, Biomacromolecules 2020, 21, 3523–3538,
  [10.1021/acs.biomac.0c00546](https://doi.org/10.1021/acs.biomac.0c00546)** — an
  ELP-*specific* temperature-dependent CG model (SOP-IDP + a semi-empirical T-dependence of
  hydrophobic interactions) that reproduces ELP cloud points vs. length and guest residue, and
  whose authors explicitly say it is "ideally suited for simulations of large-scale structures
  such as **ELP networks and hydrogels**". This is the closest published model to what this
  project is trying to do.
* **Chen & Zeng, JACS Au 2026, 6(8), 4478–4491,
  [10.1021/jacsau.6c00523](https://doi.org/10.1021/jacsau.6c00523)** — "TEA": derives
  temperature-dependent λ from hydration free energies and **bolts it onto existing one-bead
  models including CALVADOS**, validated on **16 ELP sequences** over 250–400 K with one fitted
  parameter. If a T-dependent CALVADOS is wanted without leaving the CALVADOS ecosystem, this
  is the route. *(Metadata read from the PMC record; I did not read the full paper.)*

For reference on what an ELP condensate actually is structurally: Rauscher & Pomès, eLife 2017,
6:e26526, [10.7554/eLife.26526](https://doi.org/10.7554/eLife.26526) — atomistic MD showing
the aggregated state is stabilised by the hydrophobic effect but forms **no hydrophobic core**
and stays heavily hydrated and maximally disordered. Useful as a caution against reading too
much structural order into CG condensates.

**No CALVADOS-specific ELP paper was found.** I searched PubMed and the web for CALVADOS
applied to ELPs / VPGVG phase separation 2023–2026 and found none. Recent CALVADOS extensions
are for PEG crowders ([10.1002/pro.70232](https://doi.org/10.1002/pro.70232)), phosphorylation
([10.1016/j.bpj.2025.07.001](https://doi.org/10.1016/j.bpj.2025.07.001)) and RNA
([10.1021/acs.jctc.4c01646](https://doi.org/10.1021/acs.jctc.4c01646)) — none about ELPs or
about surfaces. **Using CALVADOS2 for grafted ELP brushes is therefore an extrapolation
outside its validation set**, and should be stated as such in any write-up.

### 4.2 Ionic strength 0.19 M and pH 7.5 — yes, sensible

* Lysine side-chain pKa in folded proteins: **10.5 ± 1.1** (n = 35 measured values), Grimsley,
  Scholtz & Pace, Protein Sci. 2009, 18, 247–251,
  [10.1002/pro.19](https://doi.org/10.1002/pro.19). At pH 7.5 that is **>99.9% protonated**
  *[derived]*, so a fixed +1 charge on K is right. (CALVADOS only treats His titration
  explicitly, via Henderson–Hasselbalch with pKa 6.00 —
  [arXiv:2504.10408](https://doi.org/10.48550/arXiv.2504.10408).)
* 0.19 M gives a Debye length of 0.304/sqrt(0.19) = **0.70 nm** *[derived]*, vs 0.78 nm at
  physiological 0.15 M. Both are well inside CALVADOS's 4 nm electrostatic cutoff. Fine.
* **One caveat relevant to this project:** the chemistry is pH-sensitive in a way the model is
  blind to. Lim *et al.* 2007 got **\|G*\| 5.8 → 9.8 kPa (KV7F) and 25.8 → 45.8 kPa (KV2F)**
  going from pH 7.5 to pH 12 ([10.1021/bm061059m](https://doi.org/10.1021/bm061059m)), because
  deprotonated lysine is the nucleophile. In CALVADOS, raising the pH would change the K charge
  and *loosen* the brush — the opposite sign of effect from the experiment, where more
  crosslinks form. Do not compare pH trends between model and experiment.

### 4.3 Langevin friction 0.01 ps^-1 — dynamics are not physical

This is the CALVADOS default and is confirmed in the software paper: "Langevin integrator (by
default with time step t = 10 fs and friction coefficient γ = 0.01 ps^-1)"
([von Bülow *et al.* 2025, arXiv:2504.10408](https://doi.org/10.48550/arXiv.2504.10408)).

How far from water is that? *[derived, standard Stokes drag]* For a residue-sized bead,
radius a ≈ 0.3 nm in water (η = 1.0 mPa·s), ζ = 6πηa = 5.7×10^-12 kg/s; with m ≈ 110 Da =
1.8×10^-25 kg, the physical friction rate is ζ/m ≈ 3×10^13 s^-1 = **3×10^4 ps^-1**. The
CALVADOS default is **~10^6 times smaller**. The beads are effectively in vacuum: the model is
in the *underdamped/inertial* regime, chosen for sampling efficiency, not for kinetics. This is
the standard, well-understood trade-off of implicit-solvent CG models, but it has three
concrete consequences here:

1. **Event times in `crosslink_events.csv` are not times.** The docstring already says "the
   *rate* is not [meaningful]" — that is right, and the reason is this friction, plus the
   absence of solvent, activation barriers and the crosslinker's own diffusion.
2. **Even the *ordering* of events is only partly safe.** Ordering is robust if it is set by
   distance/topology (which pairs are near each other), and unsafe wherever it is set by how
   long a chain takes to reach somewhere, because relative diffusion rates of a short loop vs a
   whole chain are not correctly ranked by an inertial integrator without hydrodynamics.
3. **It pushes the model to the wrong Damköhler limit** — see §6.

If kinetics ever need to be interpreted, the minimal fix is to raise friction toward the
overdamped regime (γ ≳ 100 ps^-1 still lets dt = 10 fs run) and accept the sampling cost, and
to say explicitly that even then the time axis is a CG time, not laboratory time.

---

## 5. The surface

### 5.1 The current wall is essentially transparent — this is a bug-level finding

The expression is `step(z_wall-z)*0.5*(z_wall-z)^2` passed to an OpenMM `CustomExternalForce`,
whose energy is in kJ/mol with no prefactor. So **k = 1 kJ/mol/nm^2**. *[derived]*:

| penetration below z_wall | energy cost | in kT (T = 293.15 K) |
| --- | --- | --- |
| 0.1 nm | 0.005 kJ/mol | **0.002 kT** |
| 0.5 nm | 0.125 kJ/mol | 0.05 kT |
| 1.9 nm (reaching z = 0) | 1.8 kJ/mol | **0.74 kT** |

A bead can walk straight through the "surface" to z = 0 for less than one kT, and past it into
negative z for a couple of kT. **Since the box is periodic in z**, a bead that exits the bottom
reappears at the top of the box. The wall is a gentle nudge, not a solid.

### 5.2 What a solid impenetrable surface should be

Two standard choices:

**(a) Stiff one-sided harmonic** — minimal change, same functional form. To make 0.1 nm of
penetration cost > 5 kT *[derived]*:
```
0.5 * k * (0.1)^2 > 5 kT = 12.19 kJ/mol   ->   k > 2438 kJ/mol/nm^2
```
**Recommend k = 2500 kJ/mol/nm^2** (0.1 nm → 5.0 kT; 0.2 nm → 20 kT), or 5000 for a hard wall.
This is well below the CALVADOS bond constant (8033), so it will not destabilise dt = 10 fs,
but it *will* change the initial minimisation: the construction in `tools/surface.py` puts
beads on the tether plane, so check no bead starts below z_wall before raising k.

**(b) Integrated LJ 9-3 wall** — the physically principled form, obtained by integrating a
12-6 LJ fluid over a solid half-space:
```
U(z) = (3 sqrt(3) / 2) * eps * [ (sigma/z)^9 - (sigma/z)^3 ]
```
(standard result; see e.g. the SklogWiki entry on the 9-3 LJ potential, and Abraham & Singh
J. Chem. Phys. 1977 for the soft-repulsive-wall lineage — *the Abraham & Singh attribution is
from secondary sources; I did not read the original, treat as UNVERIFIED*). Truncating at the
minimum and shifting gives the purely repulsive WCA version, which is what polymer-brush CG
simulations standardly use for the grafting wall. The advantage over (a) is that it diverges,
so penetration is impossible at any temperature; the disadvantage is a new parameter (eps) and
a singularity to guard against during minimisation.

**Recommendation: (a) with k = 2500**, because it is a one-character change to an expression
that is already in every `prepare.py`, and because the brush never presses hard on the wall at
0.02 chains/nm^2 anyway. Note (b) in the write-up as the principled alternative.

Separately: the wall being purely repulsive means "how often a second lysine reaches the
surface in `free` mode is set by diffusion and by `distance`, not by any chemistry"
(`tools/surface.py` docstring) — correct, and worth keeping in the limitations list, since real
ELP coatings adsorb (ELPs are notoriously surface-active, which is the basis of the QCM-D
studies below).

### 5.3 Is 0.02 chains/nm^2 realistic?

*[derived]* 0.02 chains/nm^2 = **50 nm^2 per chain** = a square lattice spacing of
1/sqrt(0.02) = **7.07 nm** = 2×10^12 chains/cm^2.

Measured comparisons:

| system | reported | in chains/nm^2 |
| --- | --- | --- |
| silica-binding-peptide–ELP diblock brushes (B–ES40, 40× VPGSG, ~20 kDa), QCM-D + DLS + AFM; SBP saturation ≈0.5 mg/m^2, "approximate surface area per peptide on the order of **5 nm^2**"; brush height **h ≈ 14 nm** vs free-chain D_H = 9 ± 1 nm ("somewhat stretched") — [10.1021/acs.biomac.1c00067](https://doi.org/10.1021/acs.biomac.1c00067) | 5 nm^2/chain | **0.2** |
| adsorbed ELP coatings on PET, QCM-D, three sequence lengths 9.8–30.8 kDa; post-rinse **0.04–0.05 nmol/cm^2** — Srokowski & Woodhouse, J. Mater. Sci. Mater. Med. 2012, [10.1007/s10856-012-4772-6](https://doi.org/10.1007/s10856-012-4772-6) | 0.04–0.05 nmol/cm^2 | **0.24–0.30** *[derived]* |

So **0.02 chains/nm^2 is about 10× sparser than measured dense ELP layers**. That is defensible
and arguably deliberate — those measured layers are short ELPs (20–30 kDa) adsorbed or
peptide-anchored, whereas the simulated chains are 720–1040 residues (~60–90 kDa) and much
larger in footprint. Check it with the brush criterion instead of by comparing σ directly:

*[derived]* for N = 720 residues, an IDP scaling Rg ≈ 0.25 N^0.55 nm gives Rg ≈ 9.3 nm, so the
reduced grafting density Σ = σ·π·Rg² = 0.02 × π × 86 ≈ **5.4 ≫ 1**, and the graft spacing
(7.07 nm) is half the coil diameter. **The system is comfortably in the brush regime** — which
is also what the runs show (chains reaching 51–63 nm, per the `new_simulation.py` comment). So
0.02 is defensible *for these chain lengths*; it would not be for short ELPs.

The honest sentence for a report: *"0.02 chains/nm^2 (50 nm^2/chain) is ~10× below the
saturation density measured for short ELP layers by QCM-D
([10.1021/acs.biomac.1c00067](https://doi.org/10.1021/acs.biomac.1c00067),
[10.1007/s10856-012-4772-6](https://doi.org/10.1007/s10856-012-4772-6)), but for 700+ residue
chains it still corresponds to a strongly overlapping brush (Σ ≈ 5)."*

---

## 6. Time scales, and the Damköhler limit the model sits in

**Experimental gelation times** (all primary, all from §1/§3):
* THPP + ELP: **< 1 min**, "ELP hydrogels formed in less than one minute under physiological
  conditions" ([10.1021/bm061059m](https://doi.org/10.1021/bm061059m)); "< 1 min for all ELP
  block copolymers" at 200 mg/mL ([10.1021/bm7007982](https://doi.org/10.1021/bm7007982)).
* THPC + ELP: **6.7 ± 0.2 min to 27 ± 1.2 min**, tunable by ratio
  ([10.1021/bm3015279](https://doi.org/10.1021/bm3015279)).

**Chain relaxation times** are 6–9 orders of magnitude faster. For elastin-like peptides,
backbone conformational interconversion is **~1 ns** and the end-to-end vector correlation time
for poly(GVGVP) is **~500 ps at 298 K** *[these two numbers come from secondary summaries of
the elastin MD/NMR literature and I could not verify them against the primary text —
UNVERIFIED; the primary candidates are the GVG(VPGVG) MD study
[10.1529/biophysj.104.042168 / arXiv:q-bio/0401012] and Yao & Hong's solid-state NMR work]*.
Even taking a generously slow estimate of microseconds for a 700-residue chain's Rg
autocorrelation, that is still **≥10^6 times faster than gelation**.

**Therefore the physical Damköhler number is Da = τ_diffusion / τ_reaction ≪ 1: real ELP
crosslinking is strongly reaction-limited.** Chains explore many configurations between
successful bond-forming events. The model's rule — *react on first contact within a cutoff* —
is the **Da → ∞, diffusion-limited** limit. It is the opposite extreme from the experiment.

Why this matters for the project's actual question (network architecture):

* Diffusion-limited crosslinking freezes in whatever contacts happen to exist, giving a
  **more heterogeneous, more loop-rich, less equilibrated** network. Reaction-limited
  crosslinking samples configurations and gives a more uniform network at the same conversion.
* This is exactly the effect the loop literature quantifies: Wang *et al.* PRL 2016 show all
  cyclic topologies collapse onto one dimensionless formation parameter
  ([10.1103/PhysRevLett.116.188302](https://doi.org/10.1103/PhysRevLett.116.188302)), and the
  semibatch/slow-addition PNAS work shows that **slowing the effective reaction rate cuts the
  primary-loop fraction from 35% to 19% at 1 mM and raises the modulus from 220 Pa to 1300 Pa
  (a ~6× gain) for A2+B8** ([10.1073/pnas.1620985114](https://doi.org/10.1073/pnas.1620985114)).
  The kinetic regime changes the modulus by more than the chemistry does.

**Supporting simulation literature on the crossover:**
* Trautenberg, Sommer & Göritz, J. Chem. Soc. Faraday Trans. 1995, 91, 2649,
  [10.1039/FT9959102649](https://doi.org/10.1039/FT9959102649) — bond-fluctuation simulation of
  end-linked network formation showing an explicit **crossover from reaction-controlled to
  diffusion-limited kinetics as a function of the reaction rate**. This is the single most
  on-point citation for the claim.
* Duering, Kremer & Grest, J. Chem. Phys. 1994, 101, 8169,
  [10.1063/1.468202](https://doi.org/10.1063/1.468202) — the classic MD study of end-crosslinking
  kinetics and the resulting network's structure and dynamics.
* Poshusta, Bowman & Anseth, J. Biomater. Sci. Polym. Ed. 2002, 13, 797–815,
  [10.1163/156856202760197429](https://doi.org/10.1163/156856202760197429) — kinetic gelation
  simulation applied specifically to **in situ crosslinking biomaterials**, i.e. the
  application class this project is in. *(Bibliographic details verified via Semantic Scholar /
  publisher listing; abstract not read in full.)*
* Akagi, Fujimoto & Yasuda, Polym. J. 2025, 57, 1183–1194,
  [10.1038/s41428-025-01089-7](https://doi.org/10.1038/s41428-025-01089-7) — recent CG-MD
  comparison of **different crosslinking protocols** on end-linked network structure.
  *(Citation verified; full text paywalled, contents not read — UNVERIFIED.)*

**The cheap fix.** Add a per-check reaction probability `p ≪ 1` alongside the distance
criterion. `p` is the knob that moves the model from Da → ∞ toward the experimental Da ≪ 1,
and sweeping it (p = 1, 0.1, 0.01) directly tests whether the reported network architecture is
a robust result or an artefact of the diffusion-limited rule. `tools/surface.py` already has
exactly this construct (`prob` per check) for surface binding; `tools/crosslink.py` does not.
This is a small, high-value asymmetry to fix.

---

## 7. Recommended defaults, and the 4th attachment mode

### 7.1 Recommended defaults

| parameter | current | recommended | why | citation |
| --- | --- | --- | --- | --- |
| wall stiffness (`ext_force_expr` prefactor) | **1 kJ/mol/nm^2** | **2500** (5000 for a hard wall) | 0.1 nm of penetration currently costs 0.002 kT; reaching z=0 costs 0.74 kT — the surface is transparent, and the box is periodic in z so beads wrap to the top. 2500 makes 0.1 nm cost 5 kT | *[derived]*; alternative LJ 9-3 form standard for CG brush walls |
| crosslink bond k | **2000 kJ/mol/nm^2** | **250** (sweep 100 / 250 / 500) | represents 2 lysine side chains + a short spacer, contour ≈1.5 nm; FJC entropic stiffness is ~20, and k=2000 gives sigma=0.035 nm — as rigid as a peptide bond. 250 gives sigma = 0.10 nm | FJC *[derived]*; CA→NZ 0.63 nm and side-chain flexibility from [10.1002/pro.2458](https://doi.org/10.1002/pro.2458); backbone k=8033 from [arXiv:2504.10408](https://doi.org/10.48550/arXiv.2504.10408) |
| crosslink r0 | **0.6 nm** | **keep 0.6** (THPP/THPC); 0.75 glutaraldehyde, 0.9 BS3-like | sqrt(Lc·b) = 0.61 nm for a THPP-length tether — the current value is right for the right reason | *[derived]* from [10.1002/pro.2458](https://doi.org/10.1002/pro.2458) geometry |
| capture radius | **0.8 nm** | **1.2 nm** default; sweep **0.8 / 1.2 / 1.6** | hard geometric max for a THPP bridge is ≈1.6 nm (0.63+0.63+0.3); Merkley's MD-calibrated rule adds 10–25% over the hard max, and only 19% of real crosslinks satisfy the naive static criterion | [10.1002/pro.2458](https://doi.org/10.1002/pro.2458) |
| **reaction probability** | **absent (p = 1)** | **add `prob`, default 0.1; sweep 1 / 0.1 / 0.01** | react-on-first-contact is the Da→∞ diffusion-limited limit; real ELP gelation (1–27 min vs ns–µs chain relaxation) is Da ≪ 1, reaction-limited. Kinetic regime changes loop fraction 35%→19% and modulus ~6× | [10.1039/FT9959102649](https://doi.org/10.1039/FT9959102649), [10.1073/pnas.1620985114](https://doi.org/10.1073/pnas.1620985114), gel times [10.1021/bm3015279](https://doi.org/10.1021/bm3015279) |
| **junction functionality** | `valence` 1 or 2, **pairwise only** | **implicit f-functional junction, f = 3 (THPP) / 4 (THPC)** | a degree-2 node is topologically invisible; pairwise K–K bonds at any valence make paths and rings, never branch points, so they cannot build a load-bearing network | [10.1021/bm061059m](https://doi.org/10.1021/bm061059m) (f=3), [10.1021/bm3015279](https://doi.org/10.1021/bm3015279) (f=4), [10.1021/acs.macromol.5b02243](https://doi.org/10.1021/acs.macromol.5b02243) |
| Langevin friction | 0.01 ps^-1 (CALVADOS default) | **keep for structure; document loudly. Use ≥100 ps^-1 only if kinetics are ever claimed** | ~10^6× below water's Stokes friction for a residue-sized bead — the dynamics are near-inertial and the time axis is not physical | [arXiv:2504.10408](https://doi.org/10.48550/arXiv.2504.10408); *[derived]* |
| T / ionic / pH | 293.15 K / 0.19 M / 7.5 | **keep** | Lys pKa 10.5 ± 1.1 → >99.9% charged at pH 7.5; Debye length 0.70 nm, well inside the 4 nm cutoff | [10.1002/pro.19](https://doi.org/10.1002/pro.19) |
| force field | CALVADOS2 | **keep for isothermal questions; state that LCST is absent**. If the transition matters: HPS-T, the Baul ELP model, or TEA-on-CALVADOS | no T-dependence in λ | [10.1021/acscentsci.9b00102](https://doi.org/10.1021/acscentsci.9b00102), [10.1021/acs.biomac.0c00546](https://doi.org/10.1021/acs.biomac.0c00546), [10.1021/jacsau.6c00523](https://doi.org/10.1021/jacsau.6c00523) |
| grafting density | 0.02 chains/nm^2 | **keep**, report as Σ ≈ 5 not as σ | ~10× sparser than saturated short-ELP layers, but for 700+ residue chains this is still a strongly overlapping brush | [10.1021/acs.biomac.1c00067](https://doi.org/10.1021/acs.biomac.1c00067), [10.1007/s10856-012-4772-6](https://doi.org/10.1007/s10856-012-4772-6) |
| modulus reporting | *(none)* | **add the graph-based ν, µ, ξ, G_ph pipeline of §2.2**, with V = footprint × brush height | this is the project's actual deliverable, and it is a pure post-processing job on data already written | [10.1098/rspa.1976.0146](https://doi.org/10.1098/rspa.1976.0146), [10.1126/science.aag0184](https://doi.org/10.1126/science.aag0184), [arXiv:2103.16181](https://arxiv.org/abs/2103.16181) |

Target range for validation: **0.25–46 kPa** (§3).

### 7.2 What a 4th attachment/crosslinking mode would need

Two candidate designs:

**(A) Explicit crosslinker beads** — THPP/THPC molecules as diffusing particles with f reactive
arms, each arm able to bond one lysine. Requires: a new molecule species in the CALVADOS
topology (mass, radius, λ, charge — none of which exist for a phosphine, so they would have to
be invented), a bead–lysine reaction rule per arm, arm-state bookkeeping, and a concentration
parameter. Gains: (i) **crosslinker:lysine stoichiometry becomes a real variable**, which is
the one experimental axis (§3 row 5: 0.5:1 → 250 Pa, 1:1 → 2200 Pa, >1:1 → slower) that no
implicit model can touch; (ii) singly-attached "wasted" crosslinker appears naturally; (iii)
the crosslinker's own diffusion enters the kinetics. Costs: a new species with unparameterised
CG interactions in a force field that has no non-amino-acid beads except the recently added
PEG and RNA models; substantially more state to checkpoint; and a much larger search space.

**(B) Implicit f-functional junction** — when a lysine reacts, it opens a *junction* that can
accept up to f-1 further lysines within the capture radius, within a time window. No new
species: the junction is just a constraint that binds 3 (or 4) lysine beads, implemented as
f harmonic bonds to a virtual site, or as the existing pairwise machinery plus a
"these three are one junction" record. Gains: fixes the **degree-2 problem of §2.1** — which is
the actual blocker — with modest code change, reusing the existing dormant-bond/ramp
architecture. It also makes the graph analysis of §2.2 immediately meaningful, since junction
nodes finally have degree ≥3. Costs: stoichiometry stays a fitted knob, not a modelled species;
an "effective f" must be chosen rather than emerging.

**Recommendation: (B), the implicit f-functional junction, and not (A).**

Reasoning. The project's stated goal is to *qualitatively* predict network architecture and
stiffness across ELP sequences. Architecture is a topology question, and topology is decided by
node degree — which (B) fixes and (A) also fixes, but (B) fixes it without inventing CG
parameters for a phosphine that nobody has ever parameterised for CALVADOS. Every number (A)
would buy is a number the model is in no position to predict quantitatively anyway, because the
same model has no LCST (§4.1), no physical time axis (§4.3) and no solvent. Adding an
unparameterised species to a model with those three gaps buys apparent realism, not accuracy.
(B) is roughly a day of work on machinery that already exists; (A) is a research project whose
first deliverable would be "what are λ and σ for THPP?".

The upgrade path is clean: do (B) plus the reaction-probability knob from §7.1, get the modulus
pipeline of §2.2 producing numbers in the 0.25–46 kPa range, and only then ask whether the
crosslinker:lysine stoichiometry curve is worth (A).

---

## Sources, grouped

**Crosslinker chemistry / geometry**
1. Lim, Nettles, Setton & Chilkoti, *Biomacromolecules* 8:1463–1470 (2007) — THPP, f=3. [10.1021/bm061059m](https://doi.org/10.1021/bm061059m)
2. Lim, Nettles, Setton & Chilkoti, *Biomacromolecules* 9:222–230 (2008) — ELP block architectures + THPP. [10.1021/bm7007982](https://doi.org/10.1021/bm7007982)
3. Chung, Lampe & Heilshorn, *Biomacromolecules* 13:3912–3916 (2012) — THPC, f=4. [10.1021/bm3015279](https://doi.org/10.1021/bm3015279)
4. Migneault *et al.*, *BioTechniques* 37:790–802 (2004) — glutaraldehyde. [10.2144/04375RV01](https://doi.org/10.2144/04375RV01)
5. Butler, Ng & Pudney, *J. Polym. Sci. A* 41 (2003) — genipin mechanism. [10.1002/pola.10960](https://doi.org/10.1002/pola.10960)
6. Nowatzki & Tirrell, *Biomaterials* 25:1261–1267 (2004) — HMDI on aECM protein films. [10.1016/S0142-9612(03)00635-5](https://doi.org/10.1016/S0142-9612\(03\)00635-5)
7. **Merkley *et al.*, *Protein Sci.* 23:747–759 (2014) — the CA–CA distance calibration.** [10.1002/pro.2458](https://doi.org/10.1002/pro.2458)

**Network theory**
8. Flory, *Proc. R. Soc. Lond. A* 351:351–380 (1976) — phantom networks, cycle rank. [10.1098/rspa.1976.0146](https://doi.org/10.1098/rspa.1976.0146)
9. Miller & Macosko, *Macromolecules* 9:206–211 (1976). [10.1021/ma60050a004](https://doi.org/10.1021/ma60050a004)
10. Zhong, Wang, Kawamoto, Olsen & Johnson, *Science* 353:1264–1268 (2016) — RENT. [10.1126/science.aag0184](https://doi.org/10.1126/science.aag0184)
11. Kawamoto, Zhong, Wang, Olsen & Johnson, *Macromolecules* 48:8980–8988 (2015) — loops vs functionality. [10.1021/acs.macromol.5b02243](https://doi.org/10.1021/acs.macromol.5b02243)
12. Wang, Alexander-Katz, Johnson & Olsen, *PRL* 116:188302 (2016) — universal cyclic topology. [10.1103/PhysRevLett.116.188302](https://doi.org/10.1103/PhysRevLett.116.188302)
13. Zhou *et al.*, *PNAS* 109:19119–19124 (2012) — counting primary loops experimentally. [10.1073/pnas.1213169109](https://doi.org/10.1073/pnas.1213169109)
14. Semibatch loop-defect control, *PNAS* 114 (2017) — G = G0(1 − 8/3 φ), 220→1300 Pa. [10.1073/pnas.1620985114](https://doi.org/10.1073/pnas.1620985114)
15. Lin, Wang, Johnson & Olsen, *Macromolecules* 52:1685–1694 (2019). [10.1021/acs.macromol.8b01676](https://doi.org/10.1021/acs.macromol.8b01676)
16. "On the Elasticity of Polymer Model Networks Containing Finite Loops" — c_f coefficients. [arXiv:2103.16181](https://arxiv.org/abs/2103.16181)
17. Rubinstein & Panyukov, *Macromolecules* 35:6670–6686 (2002) — entanglement contribution. [10.1021/ma0203849](https://doi.org/10.1021/ma0203849)

**ELP gels (calibration)**
18. Trabbic-Carlson, Setton & Chilkoti, *Biomacromolecules* 4:572–580 (2003). [10.1021/bm025671z](https://doi.org/10.1021/bm025671z)
19. McHale, Setton & Chilkoti, *Tissue Eng.* 11:1768–1779 (2005) — transglutaminase. [10.1089/ten.2005.11.1768](https://doi.org/10.1089/ten.2005.11.1768)
20. Nettles, Chilkoti & Setton, *Adv. Drug Deliv. Rev.* 62:1479–1485 (2010) — review. [10.1016/j.addr.2010.04.002](https://doi.org/10.1016/j.addr.2010.04.002)

**Force field**
21. Tesei, Schulze, Crehuet & Lindorff-Larsen, *PNAS* 118:e2111696118 (2021). [10.1073/pnas.2111696118](https://doi.org/10.1073/pnas.2111696118)
22. **Tesei & Lindorff-Larsen, *Open Res. Europe* 2:94 (2023) — CALVADOS 2.** [10.12688/openreseurope.14967.2](https://doi.org/10.12688/openreseurope.14967.2)
23. Cao, von Bülow, Tesei & Lindorff-Larsen, *Protein Sci.* 33:e5172 (2024) — CALVADOS 3. [10.1002/pro.5172](https://doi.org/10.1002/pro.5172)
24. von Bülow *et al.* (2025) — CALVADOS software, all defaults. [arXiv:2504.10408](https://doi.org/10.48550/arXiv.2504.10408)
25. Dignon, Zheng, Kim & Mittal, *ACS Cent. Sci.* 5:821–830 (2019) — HPS-T. [10.1021/acscentsci.9b00102](https://doi.org/10.1021/acscentsci.9b00102)
26. Baul, Bley & Dzubiella, *Biomacromolecules* 21:3523–3538 (2020) — ELP-specific T-dependent CG. [10.1021/acs.biomac.0c00546](https://doi.org/10.1021/acs.biomac.0c00546)
27. Chen & Zeng, *JACS Au* 6:4478–4491 (2026) — TEA, T-dependent λ on CALVADOS, 16 ELPs. [10.1021/jacsau.6c00523](https://doi.org/10.1021/jacsau.6c00523)
28. Rauscher & Pomès, *eLife* 6:e26526 (2017) — liquid structure of elastin. [10.7554/eLife.26526](https://doi.org/10.7554/eLife.26526)
29. Grimsley, Scholtz & Pace, *Protein Sci.* 18:247–251 (2009) — Lys pKa 10.5 ± 1.1. [10.1002/pro.19](https://doi.org/10.1002/pro.19)

**Surfaces**
30. Self-assembly of ELP brushes on silica, *Biomacromolecules* 22:1966 (2021). [10.1021/acs.biomac.1c00067](https://doi.org/10.1021/acs.biomac.1c00067)
31. Srokowski & Woodhouse, *J. Mater. Sci. Mater. Med.* 24:71–84 (2012). [10.1007/s10856-012-4772-6](https://doi.org/10.1007/s10856-012-4772-6)

**Kinetics / Damköhler**
32. **Trautenberg, Sommer & Göritz, *J. Chem. Soc. Faraday Trans.* 91:2649 (1995) — reaction-controlled → diffusion-limited crossover.** [10.1039/FT9959102649](https://doi.org/10.1039/FT9959102649)
33. Duering, Kremer & Grest, *J. Chem. Phys.* 101:8169 (1994). [10.1063/1.468202](https://doi.org/10.1063/1.468202)
34. Poshusta, Bowman & Anseth, *J. Biomater. Sci. Polym. Ed.* 13:797–815 (2002). [10.1163/156856202760197429](https://doi.org/10.1163/156856202760197429)
35. Akagi, Fujimoto & Yasuda, *Polym. J.* 57:1183–1194 (2025). [10.1038/s41428-025-01089-7](https://doi.org/10.1038/s41428-025-01089-7)

### Explicitly NOT verified

* Genipin-crosslinked **ELP** modulus datapoint (§3) — no primary source read.
* Genipin adduct Nζ–Nζ span (§1.2).
* THPP P–CH2–N adduct N···N distance (§1.2) — from standard bond lengths, not a structure.
* ELP chain relaxation times ~500 ps / ~1 ns (§6) — from secondary summaries.
* Abraham & Singh 1977 as the 9-3 wall reference (§5.2) — secondary attribution.
* Akagi 2025 and Poshusta 2002 contents (citations verified, full texts not read).
* No CALVADOS–ELP paper exists as far as I could find; absence of evidence, searched
  PubMed + web for 2023–2026.
