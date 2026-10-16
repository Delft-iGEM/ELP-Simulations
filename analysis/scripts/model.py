"""The first-pass response-surface model, shared by 04_stats.py and 05_figures.py.

Full quadratic in the four coded design variables (x in [-1, +1] over the protocol
ranges): 4 linear + 4 pure quadratic + 6 two-way interaction terms + intercept = 15
parameters, fitted by OLS on log(primary response).
"""
from __future__ import annotations

import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).parent))
from common import ANALYSIS, BASE, CODINGS, RANGES, coded  # noqa: E402

PARAMS = ["concentration", "chain_length", "spacing_pct", "ionic_strength"]
SHORT = {"concentration": "C", "chain_length": "L", "spacing_pct": "S", "ionic_strength": "ion"}
LABEL = {"concentration": "grafting density (chains/nm²)", "chain_length": "chain length (repeats)",
         "spacing_pct": "lysine spacing (% of length)", "ionic_strength": "ionic strength (M)"}
TERMS = ([SHORT[p] for p in PARAMS] + [f"I({SHORT[p]}**2)" for p in PARAMS]
         + [f"{SHORT[a]}:{SHORT[b]}" for a, b in itertools.combinations(PARAMS, 2)])
RESPONSE = "inter_free_12A"


def formula(terms=TERMS, y="logy"):
    return f"{y} ~ " + " + ".join(terms)


def term_vars(term: str) -> set[str]:
    """Variables in a formula term: 'I(C**2)' -> {'C'}, 'C:ion' -> {'C', 'ion'}."""
    return set(re.findall(r"[A-Za-z]+", term.replace("I(", "("))) - {"I"}


def terms_without(param: str):
    return [t for t in TERMS if SHORT[param] not in term_vars(t)]


def load(response: str = RESPONSE, own_only: bool = False, coding: str = "linear") -> pd.DataFrame:
    """Non-excluded runs that have this response. Teammate runs carry only inter_teamA_12A,
    so they enter that fit and no other; own_only=True drops them there too. coding picks the
    design-variable scaling (CODINGS); "linear" is the frozen one."""
    m = pd.read_csv(ANALYSIS / "runs_metrics.csv")
    m = m[~m.excluded & m[response].notna()]
    if own_only:
        m = m[m.data_source == "trajectory"]
    m = m.reset_index(drop=True)
    x = CODINGS[coding](m).rename(columns=SHORT)
    df = pd.concat([m, x], axis=1)
    df["y"] = df[response]
    df["logy"] = np.log(df["y"])
    return df


def fit(df: pd.DataFrame, y="logy", terms=TERMS):
    return smf.ols(formula(terms, y), data=df).fit()


def to_coded(param: str, value, coding: str = "linear"):
    lo, hi = RANGES[param]
    v = np.asarray(value, float)
    if coding == "log":
        return (np.log(v) - np.log(lo)) / (np.log(hi) - np.log(lo)) * 2 - 1
    return (v - lo) / (hi - lo) * 2 - 1


def natural_per_coded(param: str, value, coding: str = "linear"):
    """d(natural value) / d(coded value) at a point: constant for linear, proportional to v for log."""
    lo, hi = RANGES[param]
    if coding == "log":
        return np.asarray(value, float) * (np.log(hi) - np.log(lo)) / 2
    return np.full_like(np.asarray(value, float), (hi - lo) / 2)


def predict_grid(res, px: str, py: str, fixed: dict | None = None, n=61, coding: str = "linear"):
    """Predicted response on a grid over the sampled ranges of px, py; others fixed (natural units)."""
    fixed = {**BASE, **(fixed or {})}
    xs = np.linspace(*RANGES[px], n)
    ys = np.linspace(*RANGES[py], n)
    X, Y = np.meshgrid(xs, ys)
    g = pd.DataFrame({SHORT[p]: float(to_coded(p, fixed[p], coding)) for p in PARAMS}, index=range(X.size))
    g[SHORT[px]] = to_coded(px, X.ravel(), coding)
    g[SHORT[py]] = to_coded(py, Y.ravel(), coding)
    logy = res.predict(g).to_numpy().reshape(X.shape)
    return xs, ys, logy


def gradient_coded(res, point_coded: dict, h=1e-4):
    """Numerical gradient of log(response) w.r.t. the coded variables at a point."""
    base = pd.DataFrame([point_coded])
    out = {}
    for k in point_coded:
        up, dn = base.copy(), base.copy()
        up[k] += h
        dn[k] -= h
        out[k] = float((res.predict(up).iloc[0] - res.predict(dn).iloc[0]) / (2 * h))
    return out
