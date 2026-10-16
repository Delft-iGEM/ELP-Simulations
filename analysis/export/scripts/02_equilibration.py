"""Step 2b - how much of each run to discard as equilibration (burn-in).

For every run and three observables (primary response, chain-mean Rg, 95 % mass height)
finds the burn-in t0 that maximises the effective number of independent samples
N_eff = (N - t0) / g(t0) (Chodera, JCTC 2016, 12, 1799), with g the statistical
inefficiency. Reads only the extracted time series; no summary by design parameter.

Writes analysis/qc/equilibration_per_run.csv and analysis/figures/qc_timeseries_sample.html.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import ANALYSIS, FIGURES, TIMESERIES  # noqa: E402
import style  # noqa: E402

SEED = 20261015
STEP = 20            # candidate burn-in grid (frames)


def statistical_inefficiency(x: np.ndarray) -> float:
    """g = 1 + 2 sum_t (1 - t/N) C(t)/C(0), summed until C(t) first drops <= 0 (pymbar)."""
    x = np.asarray(x, float)
    n = len(x)
    dx = x - x.mean()
    var = dx.var()
    if var == 0 or n < 4:
        return 1.0
    f = np.fft.rfft(dx, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n] / (var * np.arange(n, 0, -1))
    g = 1.0
    for t in range(1, n - 1):
        if acf[t] <= 0:
            break
        g += 2.0 * acf[t] * (1 - t / n)
    return max(1.0, g)


def detect_equilibration(x: np.ndarray, step: int = STEP):
    n = len(x)
    best = (0, 1.0, -1.0)
    for t0 in range(0, n // 2, step):
        g = statistical_inefficiency(x[t0:])
        neff = (n - t0) / g
        if neff > best[2]:
            best = (t0, g, neff)
    return best


def series(z) -> dict[str, np.ndarray]:
    defs = list(z["defs"])
    k = defs.index("B_free_lysines_only")
    c = int(np.flatnonzero(z["cutoffs_A"] == 12.0)[0])
    return {"inter_free_12A": z["counts"][:, k, c, 1].astype(float),
            "rg_mean_A": z["rg_A"].mean(1).astype(float),
            "h95_A": z["z_quantiles_A"][:, 2].astype(float)}


def main():
    inv = pd.read_csv(ANALYSIS / "run_inventory.csv")
    slugs = inv.loc[inv.status == "finished", "slug"].tolist()
    rows = []
    for s in slugs:
        z = np.load(TIMESERIES / f"{s}.npz")
        for obs, x in series(z).items():
            t0, g, neff = detect_equilibration(x)
            rows.append({"slug": s, "observable": obs, "t0_frame": t0, "t0_ns": t0 * 0.07,
                         "g_after_t0": g, "n_eff": neff, "g_full": statistical_inefficiency(x)})
    eq = pd.DataFrame(rows)
    (ANALYSIS / "qc").mkdir(exist_ok=True)
    eq.to_csv(ANALYSIS / "qc" / "equilibration_per_run.csv", index=False)
    print(eq.groupby("observable")["t0_ns"].describe(percentiles=[.5, .75, .9, .95]).round(1))
    print(eq.groupby("observable")[["g_after_t0", "n_eff"]].median().round(1))

    # Random sample of time series (fixed seed), raw + 100-frame running mean.
    rng = np.random.default_rng(SEED)
    pick = sorted(rng.choice(slugs, 12, replace=False))
    for obs, title, unit in [("inter_free_12A", "Inter-chain free K–K contacts", "contacts / frame"),
                             ("rg_mean_A", "Chain-mean radius of gyration", "Rg (Å)")]:
        fig = style.small_multiples(pick, lambda s: series(np.load(TIMESERIES / f"{s}.npz"))[obs],
                                    eq.set_index(["slug", "observable"])["t0_frame"].xs(obs, level=1),
                                    title=title, ytitle=unit)
        style.save(fig, FIGURES / f"qc_timeseries_sample_{obs}")


if __name__ == "__main__":
    main()
