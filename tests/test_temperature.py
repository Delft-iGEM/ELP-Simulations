"""The `temperature` run setting: the CSV column, the Shared fallback, the
range check, and the value that actually reaches the generated prepare.py.

The point of most of these is the default. Temperature was hardcoded at
293.15 K in the prepare.py template until this setting existed, so a CSV with
no temperature column, a blank cell, and a Shared left alone must all still
produce exactly that — otherwise every run in used-runs.csv stops being
comparable with every run after it.

Run from the worktree root:
    /home/tnartey/ELP-Simulations/.venv/bin/python -m pytest tests -q
"""

from __future__ import annotations

import re

import pytest
import typer

from tools.new_simulation import (
    DEFAULT_TEMPERATURE,
    MAX_TEMPERATURE,
    MIN_TEMPERATURE,
    PREPARE_TEMPLATE,
    check_temperature,
)
from tools.run_batch import COLUMN_ALIASES, REGISTRY_COLUMNS, Shared, _norm_header, read_runs

SEQ = "VPGVGVPGVGVPGKGVPGVGVPGVG"
HEADER = "Name;Sequence;concentration;nmolecules;steps;walltime"
ROW = f"run-a;{SEQ};0.02;4;10000;1:00:00"
SHARED = Shared()


def write(tmp_path, text: str):
    path = tmp_path / "runs.csv"
    path.write_text(text, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The default is load-bearing
# ---------------------------------------------------------------------------

def test_no_column_gives_the_historical_default(tmp_path):
    (spec,) = read_runs(write(tmp_path, f"{HEADER}\n{ROW}\n"), SHARED)
    assert spec.temperature == DEFAULT_TEMPERATURE == 293.15


def test_blank_cell_falls_back_to_shared(tmp_path):
    csv = f"{HEADER};temperature\n{ROW};\n"
    (spec,) = read_runs(write(tmp_path, csv), Shared(temperature=310.0))
    assert spec.temperature == 310.0


def test_cell_overrides_shared(tmp_path):
    csv = f"{HEADER};temperature\n{ROW};313.15\n"
    (spec,) = read_runs(write(tmp_path, csv), Shared(temperature=310.0))
    assert spec.temperature == 313.15


def test_template_default_is_unchanged_by_the_substitution():
    """The placeholder is what varies; the value written for a default run is the old literal."""
    assert "temp = __TEMPERATURE__," in PREPARE_TEMPLATE
    rendered = PREPARE_TEMPLATE.replace("__TEMPERATURE__", repr(DEFAULT_TEMPERATURE))
    assert "temp = 293.15," in rendered


# ---------------------------------------------------------------------------
# Column spelling
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("header", [
    "temperature", "Temperature", "temp", "T", "Temperature (K)",
    "temperature_K", "thermostat", "Thermostat Temperature",
])
def test_header_spellings_all_land_on_temperature(header):
    assert COLUMN_ALIASES[_norm_header(header)] == "temperature"


def test_celsius_header_is_not_silently_accepted():
    """'temperature (C)' normalises to 'temperature' — the unit note is stripped.

    So the header cannot carry the unit; the *value* check is what protects a
    Celsius number, which is test_celsius_value_is_named below.
    """
    assert _norm_header("temperature (C)") == "temperature"


# ---------------------------------------------------------------------------
# The range check
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("value", [MIN_TEMPERATURE, 293.15, 310, 400, MAX_TEMPERATURE])
def test_sensible_values_pass(value):
    assert check_temperature(value) == float(value)


@pytest.mark.parametrize("value", [MIN_TEMPERATURE - 1, MAX_TEMPERATURE + 1, 0, 1e4])
def test_values_outside_the_range_are_refused(value):
    with pytest.raises(typer.BadParameter):
        check_temperature(value)


@pytest.mark.parametrize("value,kelvin", [(20, 293.15), (37, 310.15), (4, 277.15)])
def test_celsius_value_is_named(value, kelvin):
    """A plain '20' is the one mistake worth a bespoke message, not a bounds error."""
    with pytest.raises(typer.BadParameter) as error:
        check_temperature(value)
    message = str(getattr(error.value, "message", None) or error.value)
    assert "Celsius" in message
    assert f"{kelvin:g}" in message


def test_nan_and_inf_are_refused():
    for bad in (float("nan"), float("inf"), float("-inf")):
        with pytest.raises(typer.BadParameter):
            check_temperature(bad)


def test_bad_cell_is_reported_with_its_row(tmp_path):
    csv = f"{HEADER};temperature\n{ROW};20\n"
    with pytest.raises(ValueError) as error:
        read_runs(write(tmp_path, csv), SHARED)
    assert "Celsius" in str(error.value)


def test_non_numeric_cell_is_reported(tmp_path):
    csv = f"{HEADER};temperature\n{ROW};warm\n"
    with pytest.raises(ValueError) as error:
        read_runs(write(tmp_path, csv), SHARED)
    assert re.search(r"row 2", str(error.value))


# ---------------------------------------------------------------------------
# It reaches the registry
# ---------------------------------------------------------------------------

def test_registry_has_a_temperature_column():
    assert "temperature_K" in REGISTRY_COLUMNS


def test_registry_row_records_it(tmp_path):
    from tools.run_batch import record_runs

    csv = f"{HEADER};temperature\n{ROW};313.15\n"
    (spec,) = read_runs(write(tmp_path, csv), SHARED)
    path = record_runs([(spec, "123")], csv_name="t.csv", path=tmp_path / "used.csv")
    rows = path.read_text().splitlines()
    assert "temperature_K" in rows[0]
    assert "313.15" in rows[1]


def test_old_registry_rows_get_the_historical_temperature(tmp_path):
    """A used-runs.csv written before the column existed is upgraded, not blanked.

    Every one of those runs was at 293.15 K — the template had the number
    hardcoded — so a blank would throw away something known.
    """
    from tools.run_batch import record_runs

    old_columns = [c for c in REGISTRY_COLUMNS if c != "temperature_K"]
    path = tmp_path / "used.csv"
    path.write_text(";".join(old_columns).replace(";", ",") + "\n"
                    + ",".join(["old-run", "VPGVG", "5"] + [""] * (len(old_columns) - 3)) + "\n")

    csv = f"{HEADER};temperature\n{ROW};313.15\n"
    (spec,) = read_runs(write(tmp_path, csv), SHARED)
    record_runs([(spec, "1")], path=path)

    rows = list(__import__("csv").DictReader(path.read_text().splitlines()))
    by_name = {r["name"]: r for r in rows}
    assert by_name["old-run"]["temperature_K"] == f"{DEFAULT_TEMPERATURE:g}"
    assert by_name["run-a"]["temperature_K"] == "313.15"
