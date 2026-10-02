"""CALVADOS-TEA: temperature-dependent stickiness.

Chen, F. & Zeng, X. (2026), JACS Au 6(8), 4478-4491, doi 10.1021/jacsau.6c00523.
TEA = Temperature-dependent Energetics derived from hydrAtion free energies.

Stock CALVADOS gives every residue a constant lambda (stickiness) in
``residues_CALVADOS2.csv``. TEA makes lambda a function of temperature, so that
chain collapse on heating emerges from the force field instead of having to be
put in by hand. Three chained equations:

1. Hydration free energy, Gibbs-Helmholtz integrated from T0 = 300 K (eq 2)::

       du_h(T) = (du_h(T0) - dh)*T/T0 + dh + dcp*(T*(1 - ln(T/T0)) - T0)

2. Excess pair free energy, linear in the hydration change (eq 3)::

       dG_E(T) - dG_E(T0) = kappa_i * (du_h(T) - du_h(T0))

3. Mapped onto the CG stickiness::

       lambda_i(T) = lambda_i(T0) - gamma * (dG_E(T) - dG_E(T0))

Sign
----
The paper's printed eq 6 has a PLUS with gamma positive. That contradicts the
rest of the paper: Figure 1 and the text have dG_E *decreasing* for hydrophobic
analogues as T rises, so a plus sign would lower lambda on heating and predict
UCST for ELPs, the opposite of the paper's own result. The MINUS form above is
what is implemented. It reproduces the verification targets in tests/test_tea.py
to 1e-4 for every residue whose reference lambda could be confirmed.

NOT verified against the authors' own lambda-generation code: SI archive
au6c00523_si_003.zip is not available on this machine, so the sign could not be
cross-checked against their source as intended. The agreement with the
independently computed check values is the evidence standing in for it. Worth
redoing if the archive turns up; nothing downstream is blocked on it.

Residues outside [0, 1]
-----------------------
lambda is not clamped anywhere in CALVADOS - init_ah_interactions builds a raw
OpenMM expression and passes lambda straight in - and TEA pushes several
residues out of [0, 1] inside the usable range:

  lambda > 1:  Trp above ~304 K, Tyr above ~310 K, Phe above ~326 K.
               The well deepens past the plain-LJ minimum (Phe at 370 K:
               -0.956 kJ/mol against -0.833 at lambda=1). Continuous and
               differentiable, so it is well behaved, just no longer bounded
               by LJ. None of these appear in the VPGXG sweep, but Trp and Tyr
               cross almost immediately - relevant for anything aromatic, e.g.
               tryptophan-rich AMP-ELP fusions.
  lambda < 0:  Glu below ~305 K onwards (-0.0103 at 370 K). The r > 2^(1/6)s
               branch carries a factor 4*lambda, so the attractive tail does
               not vanish, it INVERTS: Glu at 370 K is +0.0047 kJ/mol repulsive
               at 0.8 nm where it was -0.0003 attractive at 300 K. At 1.5e-3 kT
               this is thermally invisible; it is documented because it is a
               sign flip rather than a clamp.

Two different claims, often conflated
-------------------------------------
*Responsiveness* is d_lambda, how far lambda moves per kelvin. Ile leads:
+0.3130 over 300-370 K at gamma=3, against Phe's +0.2806 and Val's +0.1959.

*Absolute stickiness* is lambda(T) itself, and Phe leads at every temperature
in range (1.148 vs Ile's 0.855 at 370 K), because it starts far higher
(0.867 vs 0.542 at 300 K).

Transition temperature tracks absolute lambda, not the slope. So for the
VPGXG series: VPGFG should collapse at the LOWEST temperature, while VPGIG
shows the STEEPEST response. Both are true and they are not the same ranking.

Which positions actually respond
--------------------------------
In VPGXG, three of the five positions respond, not one. Val and Pro share the
propane analogue, so each carries d_lambda = +0.1959 over 300-370 K; only the
two glycines are near-inert at +0.0155. The guest adds its own on top.

The consequence matters for reading the sweep: every VPGXG sequence has the
same strongly responsive V-P backbone, so all four collapse on heating and
they differ only by the guest. (VPGEG)84 is NOT an inert control - Glu's own
d_lambda is -0.011, thermally invisible, but the chain it sits in still
collapses. Seeing VPGEG collapse is the backbone doing its job, not the model
failing.

What this does not touch
------------------------
sigma is left exactly as CALVADOS has it; only the lambda column changes. The
temperature-dependent dielectric is already in stock CALVADOS and is not
reimplemented here. lambda is evaluated once per run, at the run temperature,
and does not vary within a run - what the authors did (their Methods).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

T0 = 300.0  # K, the reference the tables are quoted at

# CHARMM/CGenFF validation range of the underlying hydration fits. Outside it
# the Gibbs-Helmholtz extrapolation is unconstrained by any data, so it raises
# rather than returning a number that looks usable.
T_MIN = 280.0
T_MAX = 380.0

# Outside T_MIN..T_MAX a caller may pass allow_extrapolation=True and get a
# number anyway. These are the walls that flag does NOT open: beyond them the
# eq-2 form stops being a smooth continuation of the fit and starts inventing
# curvature, because the dcp * (T(1 - ln(T/T0)) - T0) term grows without any
# data to hold it down. 240 K is below any liquid-water experiment the fits
# could have come from; 420 K is above the boiling point at 1 bar.
T_HARD_MIN = 240.0
T_HARD_MAX = 420.0

# gamma scales step 3. The paper recommends 2, but that recommendation is for
# HPS-Urry; on CALVADOS gamma=3 is better on both of their own metrics (SI
# S5/S6: single-chain MAE 0.247 vs 0.266, ELP T_cloud Pearson r 0.47 vs 0.33,
# and the gamma=2 interval spans zero). They never tested CALVADOS outside
# {2, 3}, so stay inside it.
DEFAULT_GAMMA = 3.0

# The Arg row to use. The table also carries a charmm36 Arg fit (kappa = +0.059)
# whose sign is opposite; no residue maps to it, so it is inert here and kept
# only so the data file matches the published table.
_ARG_CHARMM36 = "Arg(charmm36)"

_DATA = Path(__file__).resolve().parent


def _read(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def load_tables(hydration: Path | None = None,
                kappa: Path | None = None) -> tuple[dict, dict, dict]:
    """(hydration by group, kappa by group, one-letter residue -> group).

    Val and Pro share one entry (both modelled as propane) and Gly uses the
    Backbone/Gly entry, so the residue -> group map is many-to-one.
    """
    hyd_rows = _read(hydration or _DATA / "tea_hydration.csv")
    kap_rows = _read(kappa or _DATA / "tea_kappa.csv")

    hyd = {r["group"]: {"du_h_T0": float(r["du_h_T0_kcal_mol"]),
                        "dh": float(r["dh_kcal_mol"]),
                        "dcp": float(r["dcp_cal_mol_K"]) * 1e-3,   # cal -> kcal
                        "model_compound": r["model_compound"],
                        "fit_mae": float(r["fit_mae_kcal_mol"])}
           for r in hyd_rows}
    kap = {r["group"]: {"kappa": float(r["kappa"]), "r": float(r["r"]),
                        "fit_mae": float(r["fit_mae"])}
           for r in kap_rows}

    residue_group: dict[str, str] = {}
    for r in hyd_rows:
        for one in r["residues"].strip():
            if one in residue_group:
                raise ValueError(f"residue {one!r} is in two hydration groups")
            residue_group[one] = r["group"]

    missing = {g for g in residue_group.values()} - set(kap)
    if missing:
        raise ValueError(f"no kappa for group(s): {sorted(missing)}")
    return hyd, kap, residue_group


def check_temperature(T: float, allow_extrapolation: bool = False) -> float:
    """`T` in K, or a refusal. Outside the fitted range this is extrapolation.

    `allow_extrapolation` is the caller saying, on the record, that it wants a
    number outside the fit range and will carry the caveat itself. It is never
    the default, and it does not reach past T_HARD_MIN..T_HARD_MAX.
    """
    value = float(T)
    if not math.isfinite(value):
        raise ValueError(f"temperature must be finite, got {value!r}")
    if not T_MIN <= value <= T_MAX:
        if not allow_extrapolation:
            raise ValueError(
                f"TEA is fitted over {T_MIN:g}-{T_MAX:g} K (the CHARMM/CGenFF validation "
                f"range); {value:g} K is outside it and the Gibbs-Helmholtz form is "
                f"unconstrained there. Refusing to extrapolate. Pass "
                f"allow_extrapolation=True if that is what you mean."
            )
        if not T_HARD_MIN <= value <= T_HARD_MAX:
            raise ValueError(
                f"{value:g} K is outside {T_HARD_MIN:g}-{T_HARD_MAX:g} K, which "
                f"allow_extrapolation does not open: the eq-2 heat-capacity term is "
                f"no longer a continuation of the fit out there."
            )
    return value


def du_h(group: str | Iterable[str], T: float, hydration: dict | None = None,
         allow_extrapolation: bool = False):
    """Hydration free energy in kcal/mol at `T`, eq 2. Scalar or array."""
    check_temperature(T, allow_extrapolation)
    hyd = hydration if hydration is not None else load_tables()[0]
    groups = [group] if isinstance(group, str) else list(group)
    u0 = np.array([hyd[g]["du_h_T0"] for g in groups])
    dh = np.array([hyd[g]["dh"] for g in groups])
    dcp = np.array([hyd[g]["dcp"] for g in groups])
    out = (u0 - dh) * T / T0 + dh + dcp * (T * (1.0 - np.log(T / T0)) - T0)
    return float(out[0]) if isinstance(group, str) else out


def delta_g_excess(group: str | Iterable[str], T: float, tables=None,
                   allow_extrapolation: bool = False):
    """dG_E(T) - dG_E(T0), eq 3. Scalar or array."""
    hyd, kap, _ = tables if tables is not None else load_tables()
    groups = [group] if isinstance(group, str) else list(group)
    kappa = np.array([kap[g]["kappa"] for g in groups])
    delta = (du_h(groups, T, hyd, allow_extrapolation)
             - np.array([hyd[g]["du_h_T0"] for g in groups]))
    out = kappa * delta
    return float(out[0]) if isinstance(group, str) else out


def lambda_T(residue: str | Sequence[str], T: float, gamma: float = DEFAULT_GAMMA,
             lambda_ref: float | Sequence[float] | None = None, tables=None,
             allow_extrapolation: bool = False):
    """Stickiness at `T`, eq 6 in the minus form. Scalar in, scalar out.

    `lambda_ref` is lambda(T0) for the same residue(s) — CALVADOS's own value,
    which the caller reads from residues_CALVADOS2.csv rather than this module
    carrying a copy that could drift from it. Omitted, it is read from the
    project's residues CSV.
    """
    hyd, kap, residue_group = tables if tables is not None else load_tables()
    one_in = isinstance(residue, str) and len(residue) == 1
    residues = [residue] if one_in else list(residue)

    unknown = [r for r in residues if r not in residue_group]
    if unknown:
        raise KeyError(f"no TEA parameters for residue(s) {unknown} — "
                       f"known: {''.join(sorted(residue_group))}")

    if lambda_ref is None:
        reference = reference_lambdas()
        ref = np.array([reference[r] for r in residues])
    else:
        ref = np.atleast_1d(np.asarray(lambda_ref, dtype=float))
        if len(ref) != len(residues):
            raise ValueError(f"lambda_ref has {len(ref)} value(s) for {len(residues)} residue(s)")

    groups = [residue_group[r] for r in residues]
    out = ref - gamma * delta_g_excess(groups, T, (hyd, kap, residue_group),
                                       allow_extrapolation)
    return float(out[0]) if one_in else out


def residues_csv_path(root: Path | None = None) -> Path:
    base = Path(root) if root else _DATA.parent
    return base / "residues_CALVADOS2.csv"


def reference_lambdas(path: Path | None = None) -> dict[str, float]:
    """one-letter residue -> lambda(T0), straight out of CALVADOS's own table."""
    return {r["one"]: float(r["lambdas"]) for r in _read(path or residues_csv_path())}


def lambda_table(T: float, gamma: float = DEFAULT_GAMMA,
                 path: Path | None = None,
                 allow_extrapolation: bool = False) -> dict[str, float]:
    """Every residue in the CALVADOS table -> lambda(T).

    Residues with no TEA group (the "Z" surface anchor and the "X" tagged
    lysine this project adds) keep their reference lambda: they are tags on a
    real residue, and the tag itself has no hydration thermodynamics. Their
    underlying residue is handled on its own row.
    """
    tables = load_tables()
    reference = reference_lambdas(path)
    _, _, residue_group = tables
    out = {}
    for one, ref in reference.items():
        out[one] = (lambda_T(one, T, gamma, ref, tables, allow_extrapolation)
                    if one in residue_group else ref)
    return out


def write_residues_csv(out_path: Path, T: float, gamma: float = DEFAULT_GAMMA,
                       source: Path | None = None,
                       allow_extrapolation: bool = False) -> dict[str, float]:
    """Write a CALVADOS residues table with the lambda column at `T`.

    Everything else — sigma, MW, charge, bond length, row order — is copied
    through byte for byte. Returns the lambda vector that was written, for the
    run manifest.

    This writes a NEW file; the project's residues_CALVADOS2.csv is read and
    never modified, so a TEA run cannot change what a non-TEA run simulates.
    """
    source_path = source or residues_csv_path()
    rows = _read(source_path)
    fieldnames = list(rows[0])
    lambdas = lambda_table(T, gamma, source_path, allow_extrapolation)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            row = dict(row)
            row["lambdas"] = repr(lambdas[row["one"]])
            writer.writerow(row)
    return lambdas


def report(T: float, gamma: float = DEFAULT_GAMMA, path: Path | None = None,
           allow_extrapolation: bool = False) -> str:
    """A text block for the run log: what TEA changed and what to watch."""
    reference = reference_lambdas(path)
    new = lambda_table(T, gamma, path, allow_extrapolation)
    _, _, residue_group = load_tables()
    moved = sorted((abs(new[r] - reference[r]), r) for r in residue_group)
    lines = [f"CALVADOS-TEA: lambda re-evaluated at T = {T:g} K, gamma = {gamma:g} "
             f"(T0 = {T0:g} K); sigma untouched"]
    if not T_MIN <= float(T) <= T_MAX:
        lines.append(f"  ! EXTRAPOLATED: {float(T):g} K is outside the {T_MIN:g}-{T_MAX:g} K "
                     f"fit range. No data constrains lambda here; treat it as a "
                     f"model continuation, not a prediction.")
    biggest = ", ".join(f"{r} {reference[r]:.4f}->{new[r]:.4f}" for _, r in moved[-4:][::-1])
    lines.append(f"  largest shifts: {biggest}")
    out_of_range = {r: v for r, v in new.items() if not 0.0 <= v <= 1.0}
    if out_of_range:
        lines.append("  ! outside [0, 1]: "
                     + ", ".join(f"{r} {v:+.4f}" for r, v in sorted(out_of_range.items())))
        lines.append("    CALVADOS does not clamp lambda (interactions.py init_ah_interactions "
                     "evaluates it literally): lambda>1 pushes the plateau below the LJ "
                     "minimum, lambda<0 inverts the attractive tail.")
    return "\n".join(lines)
