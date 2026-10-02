"""Verification targets for CALVADOS-TEA (tools/tea.py).

Every number in the task's VERIFICATION TARGETS block is regenerated here.
Tolerance 1e-3 on du_h, 1e-4 on lambda, as specified.

lambda_Ile(300) was originally specified as 0.6442. That is Leu's value, picked
up from the row above it: calvados/data/residues.csv has LEU 0.6440005 directly
above ILE 0.5423624. Resolved to 0.542362, which is what both this project's
table and the one CALVADOS ships carry, and the Ile targets below follow from it.

The 370 K targets below are the corrected ones. The first set was computed with
ln(370/300) rounded near the 5th decimal, which the x370 in T*(1-ln(T/T0))
amplified into ~3e-4 on du_h, scaled per residue by dcp - hence 350 K being
clean while every residue drifted at 370 K. The corrected values reproduce to
<= 4e-6 and are pinned at 1e-5.

lambda(370) is pinned as an IDENTITY rather than against a transcribed number:
lambda(T0) read from the CSV at full precision, plus gamma*|kappa*d|. The
transcribed lambda column mixes precisions - Ile was recomputed from the full
baseline (matches to 4e-7), Phe and Glu were carried over from 4-dp baselines
(0.8672, 0.0007, which match those exactly), and Val matches neither at 0.404274
against 0.404174 from the 4-dp route, consistent with a single digit in the 4th
decimal. Pinning the identity keeps the test free of that provenance; the
transcribed values are kept below as a documented cross-check.

Run from the worktree root:
    /home/tnartey/ELP-Simulations/.venv/bin/python -m pytest tests -q
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from tools import tea

TOL_DU_H = 1e-3
TOL_LAMBDA = 1e-4
# d_lambda is a difference of two 4-dp printed targets, so it carries up to
# 1e-4 of their rounding before any arithmetic. 2e-4 is the same claim.
TOL_DLAMBDA = 2e-4
# The corrected 370 K targets are exact, so this is a floating-point guard now,
# not an allowance for arithmetic. The derived bound below is kept as a
# secondary check only.
TOL_TARGET = 1e-5
GAMMA = 3.0

TABLES = tea.load_tables()
REFERENCE = tea.reference_lambdas()

# residue group -> (du_h(T), d(du_h), kappa*d(du_h))
TARGETS_370 = {   # corrected: ln(370/300) at full precision
    "Val/Pro": (3.088668, 0.548668, -0.0652915),
    "Ile":     (3.353508, 0.613508, -0.1042964),
    "Phe":     (1.304988, 0.934988, -0.0934988),
    "Glu":    (-93.387219, 3.662781, +0.0036628),
    "Backbone/Gly": (-5.855283, 1.294717, -0.0051789),
}
TARGETS_350 = {
    "Val/Pro": (2.991583, 0.451583, -0.053738),
    "Ile":     (3.270577, 0.530577, -0.090198),
    "Phe":     (1.131664, 0.761664, -0.076166),
    "Glu":    (-94.426832, 2.623168, +0.002623),
}

# one-letter -> (lambda(300), lambda(350), lambda(370), d_lambda(300->370))
LAMBDA_TARGETS = {
    "F": (0.8672, 1.0957, 1.147696, +0.2806),
    "I": (0.542362, 0.812956, 0.855252, +0.3130),   # baseline from residues_CALVADOS2.csv
    "V": (0.2083, 0.3695, 0.404274, +0.1959),
    "E": (0.0007, -0.0072, -0.010288, -0.0110),
}


# ---------------------------------------------------------------------------
# Step 1 and 2: hydration and excess free energy
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("T,targets", [(370.0, TARGETS_370), (350.0, TARGETS_350)])
def test_du_h_and_kappa_products(T, targets):
    """The model's own outputs, pinned tight. These are exact, not rounded."""
    hyd, _, _ = TABLES
    for group, (want_u, want_d, want_kd) in targets.items():
        got_u = tea.du_h(group, T, hyd)
        got_d = got_u - hyd[group]["du_h_T0"]
        got_kd = tea.delta_g_excess(group, T, TABLES)
        assert got_u == pytest.approx(want_u, abs=TOL_TARGET), f"{group} du_h({T})"
        assert got_d == pytest.approx(want_d, abs=TOL_TARGET), f"{group} d(du_h)({T})"
        assert got_kd == pytest.approx(want_kd, abs=TOL_TARGET), f"{group} kappa*d({T})"


def test_ln_ratio_precision_is_not_rounded():
    """Guard the slip that produced the first 370 K set.

    Rounding ln(370/300) at the 5th decimal moves T*(1-ln(T/T0))-T0 by ~5e-3,
    which dcp turns into ~3e-4 on du_h. Catch it at the source.
    """
    import math
    # -7.596596, not the -7.596660 quoted alongside the corrected targets; that
    # quote is off by 6e-5 but did not reach the targets, which reproduce to 4e-6.
    exact = 370.0 * (1.0 - math.log(370.0 / 300.0)) - 300.0
    assert exact == pytest.approx(-7.596596, abs=1e-6)
    rounded = 370.0 * (1.0 - round(math.log(370.0 / 300.0), 5)) - 300.0
    assert abs(exact - rounded) > 1e-4     # the slip really is this visible


def test_du_h_at_T0_is_the_table_value():
    hyd, _, _ = TABLES
    for group, row in hyd.items():
        assert tea.du_h(group, tea.T0, hyd) == pytest.approx(row["du_h_T0"], abs=1e-12)


def test_dcp_is_converted_from_cal_to_kcal():
    """A missing 1e-3 would be the easiest way to get this wrong and still look plausible."""
    hyd, _, _ = TABLES
    rows = {r["group"]: r for r in csv.DictReader(
        open(Path(tea.__file__).parent / "tea_hydration.csv", encoding="utf-8-sig"))}
    for group, row in hyd.items():
        assert row["dcp"] == pytest.approx(float(rows[group]["dcp_cal_mol_K"]) * 1e-3)


# ---------------------------------------------------------------------------
# Step 3: lambda
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("one", sorted(LAMBDA_TARGETS))
def test_lambda_values(one):
    want_300, want_350, want_370, want_d = LAMBDA_TARGETS[one]
    ref = REFERENCE[one]
    assert ref == pytest.approx(want_300, abs=TOL_LAMBDA), (
        f"lambda({one}, 300) in residues_CALVADOS2.csv is {ref:.6f}, "
        f"the target says {want_300}"
    )
    _, kap, residue_group = TABLES
    # Secondary check: a du_h error reaches lambda multiplied by gamma*|kappa|,
    # so lambda can never be pinned tighter than the du_h tolerance allows.
    # With the corrected targets this is slack, which is the point.
    tol = max(TOL_LAMBDA, GAMMA * abs(kap[residue_group[one]]["kappa"]) * TOL_TARGET)
    for T, want in ((350.0, want_350), (370.0, want_370)):
        got = tea.lambda_T(one, T, GAMMA, ref, TABLES)
        assert got == pytest.approx(want, abs=max(tol, 5e-5)), f"lambda({one}, {T})"
    got_d = tea.lambda_T(one, 370.0, GAMMA, ref, TABLES) - ref
    assert got_d == pytest.approx(want_d, abs=TOL_DLAMBDA)


@pytest.mark.parametrize("one", sorted("FIVE"))
def test_lambda_is_the_baseline_plus_the_model_term(one):
    """lambda(370) pinned as an identity, free of any transcribed value's rounding.

    lambda(T0) comes from the CSV at full precision and kappa*d is pinned to 1e-5
    above, so this is the only form of the lambda target that is exactly
    determined.
    """
    group = TABLES[2][one]
    want_kd = TARGETS_370[group][2]
    expected = REFERENCE[one] - GAMMA * want_kd
    # kappa*d is pinned to TOL_TARGET, and lambda multiplies it by gamma, so the
    # identity holds exactly to gamma*TOL_TARGET and no tighter.
    assert tea.lambda_T(one, 370.0, GAMMA, REFERENCE[one], TABLES) == pytest.approx(
        expected, abs=GAMMA * TOL_TARGET)


@pytest.mark.parametrize("one", sorted(LAMBDA_TARGETS))   # Ile included: its SHIFT is fine
def test_lambda_shift_is_independent_of_the_reference(one):
    """d_lambda is what TEA contributes; it must not depend on lambda(T0) at all.

    This is the part of the Ile target that DOES reproduce, which is how we know
    the conflict is in the baseline and not in the model.
    """
    want_d = LAMBDA_TARGETS[one][3]
    for ref in (0.0, 0.5, REFERENCE[one]):
        got = tea.lambda_T(one, 370.0, GAMMA, ref, TABLES) - ref
        assert got == pytest.approx(want_d, abs=TOL_DLAMBDA)


def test_ile_baseline_is_not_leu():
    """Regression: 0.6442 is Leu's lambda, one row above Ile in the residues table.

    The two sit adjacent in calvados/data/residues.csv (LEU 0.6440005, then ILE
    0.5423624) and share MW 113.16 and sigma 0.618, so a row slip between them
    is easy to make and hard to see. Pin both.
    """
    assert REFERENCE["I"] == pytest.approx(0.542362, abs=1e-6)
    assert REFERENCE["L"] == pytest.approx(0.644001, abs=1e-6)
    assert REFERENCE["I"] != pytest.approx(REFERENCE["L"], abs=1e-3)


@pytest.mark.parametrize("T", [310.0, 350.0, 370.0])
def test_delta_lambda_is_baseline_independent_for_every_residue(T):
    """TEA's contribution cannot depend on lambda(T0) — for any residue, any baseline.

    This is what made the Ile dispute decidable: the ranking is set by d_lambda,
    so a wrong baseline could not have moved it. Pinned so the next baseline
    argument cannot reach the ranking either.
    """
    _, _, residue_group = TABLES
    for one in sorted(residue_group):
        shifts = [tea.lambda_T(one, T, GAMMA, ref, TABLES) - ref
                  for ref in (-1.0, 0.0, 0.25, REFERENCE[one], 1.0, 7.5)]
        assert max(shifts) - min(shifts) < 1e-12, f"{one} at {T} K depends on its baseline"


def test_ranking_is_I_then_F_then_V():
    """Guest-residue contribution to collapse: I > F > V >> E ~ 0."""
    shift = {one: tea.lambda_T(one, 370.0, GAMMA, 0.0, TABLES) for one in "IFVE"}
    assert shift["I"] > shift["F"] > shift["V"] > abs(shift["E"])
    assert shift["E"] == pytest.approx(-0.011, abs=2e-3)


def test_glu_signal_is_below_its_own_fit_error():
    """E is flat, not anti-responsive: the signal is smaller than the fits it passes through."""
    hyd, kap, _ = TABLES
    kappa_signal = abs(tea.delta_g_excess("Glu", 370.0, TABLES))
    assert kappa_signal == pytest.approx(0.003663, abs=TOL_DU_H)
    assert kappa_signal < hyd["Glu"]["fit_mae"]
    assert kappa_signal < kap["Glu"]["fit_mae"]


# ---------------------------------------------------------------------------
# Range, tables, and the written file
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("T", [tea.T_MIN, 300.0, 350.0, tea.T_MAX])
def test_in_range_temperatures_are_accepted(T):
    assert tea.check_temperature(T) == T


@pytest.mark.parametrize("T", [279.9, 380.1, 0.0, 1000.0])
def test_out_of_range_temperatures_raise(T):
    with pytest.raises(ValueError, match="Refusing to extrapolate"):
        tea.check_temperature(T)
    with pytest.raises(ValueError):
        tea.du_h("Phe", T)


def test_val_and_pro_share_an_entry_and_gly_uses_backbone():
    _, _, residue_group = TABLES
    assert residue_group["V"] == residue_group["P"] == "Val/Pro"
    assert residue_group["G"] == "Backbone/Gly"


def test_cgenff_arg_is_used_not_charmm36():
    _, kap, residue_group = TABLES
    assert kap[residue_group["R"]]["kappa"] == pytest.approx(-0.014)
    assert "Arg(charmm36)" in kap          # present for reference
    assert "Arg(charmm36)" not in residue_group.values()   # but never selected


def test_every_standard_residue_has_parameters():
    _, _, residue_group = TABLES
    assert set("ACDEFGHIKLMNPQRSTVWY") <= set(residue_group)


def test_vectorised_matches_scalar():
    residues = list("ACDEFGHIKLMNPQRSTVWY")
    refs = np.array([REFERENCE[r] for r in residues])
    vector = tea.lambda_T(residues, 340.0, GAMMA, refs, TABLES)
    scalar = [tea.lambda_T(r, 340.0, GAMMA, REFERENCE[r], TABLES) for r in residues]
    assert np.allclose(vector, scalar, atol=1e-12)


def test_at_T0_lambda_is_exactly_the_reference():
    """The identity that makes TEA safe to leave on at 300 K."""
    table = tea.lambda_table(tea.T0, GAMMA)
    for one, ref in REFERENCE.items():
        assert table[one] == pytest.approx(ref, abs=1e-12)


def test_written_csv_changes_only_lambda(tmp_path):
    out = tmp_path / "residues_TEA.csv"
    tea.write_residues_csv(out, 350.0, GAMMA)
    before = list(csv.DictReader(open(tea.residues_csv_path(), encoding="utf-8-sig")))
    after = list(csv.DictReader(open(out, encoding="utf-8-sig")))
    assert [r["one"] for r in before] == [r["one"] for r in after]
    for b, a in zip(before, after):
        for column in b:
            if column == "lambdas":
                continue
            assert a[column] == b[column], f"{b['one']}: {column} changed"


def test_written_csv_is_loadable_by_calvados(tmp_path):
    import pandas as pd
    out = tmp_path / "residues_TEA.csv"
    tea.write_residues_csv(out, 370.0, GAMMA)
    frame = pd.read_csv(out).set_index("one")
    assert frame.loc["F", "lambdas"] > 1.0      # the >1 case really does reach the file
    assert frame.loc["E", "lambdas"] < 0.0      # and the <0 case
    assert frame["sigmas"].notna().all()


def test_source_residues_csv_is_never_modified(tmp_path):
    before = tea.residues_csv_path().read_bytes()
    tea.write_residues_csv(tmp_path / "out.csv", 370.0, GAMMA)
    assert tea.residues_csv_path().read_bytes() == before


def test_anchor_and_tagged_residues_keep_their_reference_lambda():
    """"Z" and "X" are tags this project adds; they have no hydration analogue."""
    _, _, residue_group = TABLES
    table = tea.lambda_table(370.0, GAMMA)
    for one in REFERENCE:
        if one not in residue_group:
            assert table[one] == pytest.approx(REFERENCE[one], abs=1e-12)


# ---------------------------------------------------------------------------
# The two rankings, which are not the same ranking
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("T", [300.0, 330.0, 350.0, 370.0])
def test_absolute_stickiness_ranks_F_above_I_above_V(T):
    """Transition temperature tracks lambda(T), where Phe leads at every T.

    The other ranking — responsiveness, d_lambda, where Ile leads — is
    test_ranking_is_I_then_F_then_V. Conflating them is the easy mistake: Phe
    should collapse at the lowest temperature, Ile should respond most steeply.
    """
    table = tea.lambda_table(T, GAMMA)
    assert table["F"] > table["I"] > table["V"] > table["E"]


def test_phe_leads_absolute_but_ile_leads_response():
    """Both at once, stated as one assertion so the distinction cannot drift."""
    hot = tea.lambda_table(370.0, GAMMA)
    shift = {one: hot[one] - REFERENCE[one] for one in "FIV"}
    assert hot["F"] > hot["I"]            # absolute: Phe
    assert shift["I"] > shift["F"]        # response: Ile


def test_val_and_pro_respond_equally_so_the_backbone_is_not_inert():
    """Three of five VPGXG positions respond; only the two glycines are near-inert.

    Val and Pro share the propane analogue, so every VPGXG sequence carries the
    same strongly responsive backbone. (VPGEG)84 therefore collapses too — it is
    not an inert control.
    """
    hot = tea.lambda_table(370.0, GAMMA)
    d_val = hot["V"] - REFERENCE["V"]
    d_pro = hot["P"] - REFERENCE["P"]
    d_gly = hot["G"] - REFERENCE["G"]
    assert d_val == pytest.approx(d_pro, abs=1e-12)
    assert d_val == pytest.approx(0.1959, abs=TOL_DLAMBDA)
    assert d_gly == pytest.approx(0.0155, abs=TOL_DLAMBDA)
    assert d_val > 10 * d_gly             # backbone response is not a rounding error


# ---------------------------------------------------------------------------
# The generated prepare.py must be valid Python
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("template_name", ["new_simulation", "tea_sweep", "tea_pairs"])
def test_generated_prepare_py_parses(template_name, tmp_path):
    """Regression: writing tea.yaml from a template needs doubled backslashes.

    The template is itself a Python string, so a "\\n" in its source becomes a
    real newline in the generated file — producing an unterminated string
    literal that only shows up at prepare time, after the run folder is written.
    It bit twice (tools/tea_sweep.py, then tools/new_simulation.py), so parse
    every template's output here instead.
    """
    import ast
    import importlib

    module = importlib.import_module(f"tools.{template_name}")
    template = getattr(module, "PREPARE_TEMPLATE", None) or getattr(module, "PREPARE", None)
    if template is None:
        template = getattr(module, "FREE_PREPARE")

    if "__TEMPERATURE__" in template:        # new_simulation's placeholder style
        rendered = template
        for placeholder, value in (
            ("__MODE_NOTE__", ""), ("__MODE__", "brush"), ("__SURFACE__", "None"),
            ("__SPACING_NOTE__", ""), ("__SPACING__", "7.0"), ("__MARGIN_NOTE__", ""),
            ("__MARGIN__", "0.0"), ("__CROSSLINK_NOTE__", ""), ("__CROSSLINK__", "None"),
            ("__NX__", "4"), ("__NY__", "4"), ("__L_X__", "28.0"), ("__L_Y__", "28.0"),
            ("__Z_HEIGHT__", "50.0"), ("__N_SAVE__", "7000"), ("__N_FRAMES__", "100"),
            ("__Z_WALL__", "1.9"), ("__WALL_K__", "5000.0"), ("__Z_ANCHOR__", "2.5"),
            ("__SEQ_NAME__", "x"), ("__SEQUENCE__", "VPGVG"), ("__NMOL__", "4"),
            ("__PLATFORM__", "CPU"), ("__PARTITION__", "gpu-a100"),
            ("__WALLTIME__", "1:00:00"), ("__CPU_PER_TASK__", "18"),
            ("__TEMPERATURE__", "293.15"), ("__TEA__", "{'gamma': 3.0}"),
        ):
            rendered = rendered.replace(placeholder, value)
    else:                                     # str.format style
        rendered = template.format(
            label="x", source="y", sequence="VPGVG", temperature=350.0, gamma=3.0,
            box=60.0, n_save=7000, n_frames=100, seed=1, ionic=0.19, ph=7.5,
            cutoff_lj=2.0, cutoff_yu=4.0, friction=0.01,
            partition="gpu-a100", walltime="1:00:00", cpus="18")

    ast.parse(rendered)                      # the assertion
    assert "tea.yaml" in rendered
    assert '\\n' in rendered or "\\n" in repr(rendered)
