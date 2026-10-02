"""The 270 K runs extrapolate below TEA's fit range. That must be opt-in only."""

from __future__ import annotations

import re
from collections import Counter

import pytest

from tools import tea
from tools.tea_sequences import all_terminal, original, verify, vpgig_only


def test_default_still_refuses_below_the_fit_range():
    """No ordinary call can wander outside 280-380 K by accident."""
    for T in (270.0, 279.9, 380.1, 400.0):
        with pytest.raises(ValueError, match="Refusing to extrapolate"):
            tea.check_temperature(T)
        with pytest.raises(ValueError, match="Refusing to extrapolate"):
            tea.lambda_T("I", T)


def test_inside_the_fit_range_the_flag_changes_nothing():
    for T in (280.0, 300.0, 350.0, 380.0):
        for r in "FIVGEKRD":
            assert tea.lambda_T(r, T) == tea.lambda_T(r, T, allow_extrapolation=True)


def test_the_flag_opens_270_k():
    value = tea.lambda_T("I", 270.0, allow_extrapolation=True)
    assert 0.0 < value < tea.lambda_T("I", 280.0)


def test_hard_walls_are_not_opened_by_the_flag():
    for T in (239.9, 420.1, 1000.0):
        with pytest.raises(ValueError, match="allow_extrapolation does not open"):
            tea.check_temperature(T, allow_extrapolation=True)


@pytest.mark.parametrize("T", [240.0, 270.0, 420.0])
def test_extrapolated_lambda_stays_finite_and_monotone(T):
    """The continuation must not blow up or fold back on itself."""
    for r in "FIVPGEKRD":
        value = tea.lambda_T(r, T, allow_extrapolation=True)
        assert abs(value) < 5.0, (r, T, value)
    # lambda rises with T for every hydrophobe across the whole opened range
    for r in "FIVP":
        series = [tea.lambda_T(r, t, allow_extrapolation=True)
                  for t in (240.0, 270.0, 300.0, 350.0, 420.0)]
        assert series == sorted(series), (r, series)


def test_report_announces_the_extrapolation():
    text = tea.report(270.0, 3.0, allow_extrapolation=True)
    assert "EXTRAPOLATED" in text
    assert "not a prediction" in text
    assert "EXTRAPOLATED" not in tea.report(300.0, 3.0)


def test_written_table_at_270_is_a_real_residues_csv(tmp_path):
    out = tmp_path / "residues_TEA.csv"
    lambdas = tea.write_residues_csv(out, 270.0, 3.0, allow_extrapolation=True)
    import csv
    rows = list(csv.DictReader(open(out)))
    source = list(csv.DictReader(open(tea.residues_csv_path())))
    assert len(rows) == len(source)
    assert [r["one"] for r in rows] == [r["one"] for r in source]
    # sigma and charge copied through untouched; only lambdas moved
    for new, old in zip(rows, source):
        assert new["sigmas"] == old["sigmas"]
        assert float(new["lambdas"]) == pytest.approx(lambdas[new["one"]])


def test_write_at_270_without_the_flag_refuses(tmp_path):
    with pytest.raises(ValueError, match="Refusing to extrapolate"):
        tea.write_residues_csv(tmp_path / "x.csv", 270.0, 3.0)


# --- the designed sequences -------------------------------------------------

def test_all_terminal_is_a_permutation_of_the_original():
    orig, term = original(), all_terminal()
    assert Counter(orig) == Counter(term)
    assert len(orig) == len(term) == 357
    assert orig != term


def test_all_terminal_puts_every_rgd_and_lysine_near_an_end():
    term = all_terminal()
    n = len(term)
    rgd = [m.start() for m in re.finditer("RGD", term)]
    lys = [i for i, c in enumerate(term) if c == "K"]
    assert len(rgd) == 4 and len(lys) == 4
    # every motif sits in the outer fifth of the chain at one end or the other
    for position in rgd + lys:
        assert position < n / 5 or position > 4 * n / 5, (position, n)
    # two at each end, not all four at one
    assert sum(p < n / 2 for p in rgd) == 2
    assert sum(p < n / 2 for p in lys) == 2


def test_the_original_keeps_its_motifs_in_the_interior():
    """Guards the contrast: if the original were already terminal there is no test."""
    orig = original()
    interior = [m.start() for m in re.finditer("RGD", orig)
                if len(orig) / 5 <= m.start() <= 4 * len(orig) / 5]
    assert interior, "the original should have RGD away from both ends"


def test_vpgig_only_is_bare_and_the_same_length():
    bare = vpgig_only()
    assert len(bare) == len(original()) == 357
    assert "K" not in bare and "RGD" not in bare
    assert set(bare) == set("VPGI")
    assert bare.startswith("VPGIGVPGIG") and bare.endswith("VPGIGVP")


def test_verify_passes():
    out = verify()
    assert out["composition_identical"] and out["sequences_differ"]
    assert out["vpgig_only_has_no_K_or_RGD"]
