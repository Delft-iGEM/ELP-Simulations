"""Thermoresponse of the 357-aa RGD construct: the g6-tsw temperature sweep.

    .venv/bin/python tools/tea_tsweep_analysis.py               # analyse + plot
    .venv/bin/python tools/tea_tsweep_analysis.py --plot-only   # re-plot from the CSVs

Runs: g6-tsw-275k .. 380k (tools/tea_tsweep.py, same seed at every T; 275 K is
TEA extrapolated below its 280 K fit range) and g6-orig-285k / -370k (same
sequence and settings, another seed; open squares, a seed-to-seed check).
Unfinished runs are skipped. 400 ns, TEA gamma 6, the first 50 ns discarded.

Per run: mean Rg with 20-block SEM, RGD SASA per site (tools/tea_g6_analysis),
the fraction of frames with Rg < Rg_theta, the Rg histogram, and the collapsed
dwells / tau of tools/tea_lysramp_dwells (hysteresis 0.95 / 1.05 Rg_theta,
minimum 3 frames, censored geometric tau, 1000-resample bootstrap CI).

Transition temperatures are read off by linear interpolation between the two
sampled temperatures that bracket the crossing: T_theta where mean Rg crosses
Rg_theta, T_half where the collapsed fraction crosses 0.5. With 5 K spacing
and one run per temperature they are estimates to within a few K, not more.

Writes, in figures/g6-tsw/: g6-tsw.csv, g6-tsw-dwells.csv,
g6-tsw-dwell-summary.csv, g6-tsw-rg-hist.csv, g6-tsw.json and
g6-tsw-{rg,collapse,tau,dwells-strip,rgd-sasa,rg-distribution}.png.
"""

from __future__ import annotations

import argparse
import csv
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_lysramp_dwells import (COLD, DISCARD_NS, GRID, HIGH, HOT, INK, LOW,
                                      MIN_FRAMES, MUTED, N_BOOT, SEED, SMALL_N, SURFACE,
                                      TAU_FLOOR, analyse, rg_series, summarise)
from tools.tea_tsweep import TEMPERATURES

SWEEP = [(int(T), f"g6-tsw-{int(T)}k", "sweep") for T in TEMPERATURES]
ENDS = [(285, "g6-orig-285k", "end"), (370, "g6-orig-370k", "end")]
OUT_DIR = Path("figures/g6-tsw")
BINS = np.linspace(2.5, 9.5, 71)                 # nm, Rg histogram edges
LINE = "#5a5f66"                                 # one series: neutral ink, T on x


def both(name: str) -> dict | None:
    from tools.tea_g6_analysis import analyse_run

    row = analyse_run(name)
    if row is None:
        return None
    rg, _, _ = rg_series(name)
    hist = np.histogram(rg, bins=BINS, density=True)[0]
    return {"rg": row, "dw": analyse(name), "hist": hist.tolist()}


def crossing(T: np.ndarray, y: np.ndarray, level: float) -> float:
    """First T where y crosses `level`, linear between bracketing samples."""
    for i in range(len(T) - 1):
        a, b = y[i] - level, y[i + 1] - level
        if a == 0:
            return float(T[i])
        if a * b < 0:
            return float(T[i] + (T[i + 1] - T[i]) * a / (a - b))
    return float("nan")


# ---- figures ---------------------------------------------------------------------

def _axes(plt, height: float = 4.6, width: float = 7.4):
    fig, ax = plt.subplots(figsize=(width, height))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("Temperature (K)", fontsize=9.5, color=INK)
    ax.set_xlim(271, 384)
    ax.set_xticks(range(275, 381, 10))
    return fig, ax


def _finish(fig, plt, ax, title: str, caption: str, out: Path, extra=(),
            line: bool = True, slide: bool = False) -> None:
    """slide: transparent, no title, no caption, sweep only (no other-seed points)."""
    if slide:
        handles = [plt.Line2D([], [], color=LINE, lw=2 if line else 0, marker="o", ms=6.5,
                              mec=SURFACE, label="RGD + sequence"),
                   *extra]
    else:
        handles = [plt.Line2D([], [], color=LINE, lw=2 if line else 0, marker="o", ms=6.5,
                              mec=SURFACE, label="sweep, one seed"),
                   plt.Line2D([], [], ls="", marker="s", ms=6.5, mfc=SURFACE, mec=LINE, mew=1.6,
                              label="g6-orig, other seed"), *extra]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=len(handles), frameon=False, fontsize=8.5, labelcolor=INK)
    if slide:
        fig.tight_layout()
        fig.savefig(out, dpi=200, transparent=True)
        plt.close(fig)
        print(f"wrote {out}")
        return
    ax.set_title(title, loc="left", fontsize=10.5, color=INK, fontweight="bold")
    n_lines = caption.count("\n") + 1
    fig.text(0.01, 0.01, caption, fontsize=7.8, color=MUTED, ha="left", va="bottom",
             linespacing=1.5)
    h = fig.get_size_inches()[1]
    fig.tight_layout(rect=(0, (0.2 * n_lines + 0.08) / h, 1, 1))
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def _series(ax, rows, key, sem=None, ends: bool = True):
    sw = sorted((r for r in rows if r["kind"] == "sweep"), key=lambda r: r["T"])
    T = np.array([r["T"] for r in sw])
    y = np.array([r[key] for r in sw])
    ax.plot(T, y, color=LINE, lw=2, zorder=2)
    e = np.array([r[sem] for r in sw]) if sem else None
    ax.errorbar(T, y, yerr=e, fmt="o", ms=6.5, color=LINE, elinewidth=1, capsize=0,
                mec=SURFACE, mew=1.2, zorder=3)
    for r in rows if ends else ():
        if r["kind"] == "end":
            ax.errorbar(r["T"], r[key], yerr=r[sem] if sem else None, fmt="s", ms=6.5,
                        color=LINE, mfc=SURFACE, mew=1.6, elinewidth=1, capsize=0, zorder=3)


def _vline(ax, T: float, label: str):
    if np.isfinite(T):
        ax.axvline(T, color=HOT, lw=1, ls=(0, (4, 2)), zorder=1)
        ax.text(T + 0.6, 0.97, label, transform=ax.get_xaxis_transform(), ha="left",
                va="top", fontsize=8, color=INK)


def plot(rows, events, summary, hist, m, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    out_dir.mkdir(parents=True, exist_ok=True)
    theta = m["rg_theta_nm"]
    window = f"{m['analysis_window_ns'][0]:g}–{m['analysis_window_ns'][1]:g} ns"
    rule = (f"Collapsed: Rg < Rg_θ = 0.21·N^0.5 = {theta:.3f} nm (N = {m['n_residues']}); "
            f"dwells with hysteresis (start below {LOW:g}·Rg_θ, end above {HIGH:g}·Rg_θ),\n"
            f"minimum {MIN_FRAMES} frames ({m['min_dwell_ns']:.2f} ns); {window} analysed.")
    one = "One run per temperature."

    # 1. Rg
    fig, ax = _axes(plt)
    _series(ax, rows, "rg_nm", "rg_sem_nm", ends=False)
    ax.axhline(theta, color=MUTED, lw=1, ls=(0, (1, 2)), zorder=1)
    ax.text(0.01, theta, f"θ-state Rg {theta:.2f} nm", transform=ax.get_yaxis_transform(),
            va="bottom", ha="left", fontsize=8, color=MUTED)
    _vline(ax, m["T_theta_K"], f"Rg = Rg_θ at ≈ {m['T_theta_K']:.0f} K")
    ax.set_ylabel("Radius of gyration (nm)", fontsize=9.5, color=INK)
    _finish(fig, plt, ax, "Radius of gyration vs temperature, 357-aa RGD construct",
            f"Mean Rg over {window}; bar = SEM from 20 block averages. {one}\n"
            "Crossing temperature: linear interpolation between the bracketing sweep points.",
            out_dir / "g6-tsw-rg.png", slide=True)

    # 2. collapse fraction
    fig, ax = _axes(plt)
    _series(ax, rows, "collapse_fraction")
    ax.axhline(0.5, color=MUTED, lw=1, ls=(0, (1, 2)), zorder=1)
    _vline(ax, m["T_half_K"], f"half the time collapsed at ≈ {m['T_half_K']:.0f} K")
    ax.set_ylim(-0.03, 1.0)
    ax.set_ylabel("Fraction of time collapsed", fontsize=9.5, color=INK)
    _finish(fig, plt, ax, "Fraction of time collapsed vs temperature",
            f"Fraction of the {window} frames with Rg < Rg_θ = {theta:.3f} nm (no hysteresis, "
            f"no minimum length). {one}", out_dir / "g6-tsw-collapse.png")

    # 3. tau, points only
    fig, ax = _axes(plt)
    dt = m["frame_ns"]
    for r in rows:
        if not np.isfinite(r["tau_ns"]):
            continue
        small = r["n_dwells"] < SMALL_N
        mk = "s" if r["kind"] == "end" else "o"
        ax.plot([r["T"]] * 2, [max(r["tau_lo_ns"], TAU_FLOOR), max(r["tau_hi_ns"], TAU_FLOOR)],
                color=LINE, lw=1.3, zorder=2, solid_capstyle="butt")
        y = max(r["tau_ns"], TAU_FLOOR)
        ax.plot(r["T"], y, "v" if r["tau_ns"] < TAU_FLOOR else mk, ms=6.5, zorder=3,
                clip_on=False, mfc=SURFACE if (small or r["kind"] == "end") else LINE,
                mec=LINE, mew=1.5)
    ax.set_yscale("log")
    ax.set_ylim(TAU_FLOOR, 3.0)
    ax.axhspan(TAU_FLOOR, dt, color=GRID, alpha=0.6, lw=0, zorder=0)
    ax.text(0.01, TAU_FLOOR * 1.25, f" below one frame ({dt:.2f} ns): not resolved",
            transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8,
            color=MUTED)
    ax.set_ylabel("Dwell time τ (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, ax, "Collapsed-state dwell time τ vs temperature",
            rule + f"\nτ: censored geometric ML; bar = 95 % bootstrap CI ({N_BOOT} resamples); "
                   f"open circle = fewer than {SMALL_N} dwells; ▼ = τ of 0.",
            out_dir / "g6-tsw-tau.png",
            [plt.Line2D([], [], ls="", marker="o", ms=6.5, mfc=SURFACE, mec=MUTED, mew=1.5,
                        label=f"n < {SMALL_N} dwells")], line=False)

    # 4. every dwell, median + IQR
    rng = np.random.default_rng(SEED)
    fig, ax = _axes(plt, 4.8)
    by = {s["run"]: s for s in summary}
    for r in rows:
        s = by[r["run"]]
        d = np.array([e["duration_ns"] for e in events
                      if e["run"] == r["run"] and e["rule"] == "hyst"])
        if len(d):
            ax.scatter(r["T"] + rng.uniform(-1.2, 1.2, len(d)), d, s=5, color=LINE,
                       alpha=0.25, linewidths=0, zorder=2)
        if s["n_dwells"]:
            small = s["n_dwells"] < SMALL_N
            ax.plot([r["T"]] * 2, [s["q25_ns"], s["q75_ns"]], color=INK, lw=2, zorder=3,
                    solid_capstyle="butt")
            ax.plot(r["T"], s["median_ns"], "s" if r["kind"] == "end" else "o", ms=6.5,
                    zorder=4, mfc=SURFACE if small else INK, mec=INK, mew=1.6)
    ax.set_yscale("log")
    ax.set_ylabel("Collapsed dwell (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, ax, "Collapsed-state dwells vs temperature",
            rule + "\nOne faint point per dwell; marker = median, bar = IQR; open = fewer "
                   f"than {SMALL_N} dwells.", out_dir / "g6-tsw-dwells-strip.png", line=False)

    # 5. RGD SASA
    fig, ax = _axes(plt)
    _series(ax, rows, "rgd_sasa_per_site_nm2", "rgd_sasa_sem_per_site_nm2")
    ax.set_ylabel("SASA per RGD (nm²)", fontsize=9.5, color=INK)
    _finish(fig, plt, ax, "RGD solvent exposure vs temperature",
            f"SASA of the R, G and D beads, summed and divided by the 4 sites; every 10th "
            f"frame of {window};\nbar = 20-block SEM. {one}", out_dir / "g6-tsw-rgd-sasa.png")

    # 6. Rg distribution per temperature, sequential single hue
    fig, ax = _axes(plt, 4.8, 7.8)
    cmap = LinearSegmentedColormap.from_list("seq", [SURFACE, "#9db4e6", COLD, "#1d3570"])
    temps = sorted(hist, key=int)
    # one 5 K column per temperature; temperatures not sampled stay blank
    grid = np.arange(int(min(temps, key=int)), int(max(temps, key=int)) + 1, 5)
    edges_T = np.append(grid - 2.5, grid[-1] + 2.5)
    Z = np.full((len(BINS) - 1, len(grid)), np.nan)
    for t in temps:
        Z[:, list(grid).index(int(t))] = hist[t]
    Z = np.ma.masked_invalid(Z)
    pc = ax.pcolormesh(edges_T, BINS, Z, cmap=cmap, shading="flat", zorder=1)
    ax.grid(False)
    ax.axhline(theta, color=HOT, lw=1.2, ls=(0, (4, 2)), zorder=2)
    ax.text(0.01, theta, f" Rg_θ {theta:.2f} nm", transform=ax.get_yaxis_transform(),
            va="bottom", ha="left", fontsize=8, color=INK)
    cb = fig.colorbar(pc, ax=ax, pad=0.015)
    cb.set_label("Probability density (1/nm)", fontsize=8.5, color=INK)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_edgecolor(GRID)
    ax.set_xlim(edges_T[0], edges_T[-1])
    ax.set_ylim(2.8, 9.0)
    ax.set_ylabel("Radius of gyration (nm)", fontsize=9.5, color=INK)
    ax.set_title("Rg distribution at each temperature", loc="left", fontsize=10.5,
                 color=INK, fontweight="bold")
    fig.text(0.01, 0.01, f"Histogram of Rg over {window}, every frame, one column per "
             "temperature; the one-seed sweep only. Below 280 K TEA is extrapolated.",
             fontsize=7.8, color=MUTED)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = out_dir / "g6-tsw-rg-distribution.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


# ---- main ------------------------------------------------------------------------

def _write(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def _read(path: Path, text=("run", "kind", "rule", "tea_gamma")) -> list[dict]:
    out = []
    for r in csv.DictReader(open(path)):
        r = {k: (v if k in text else float(v)) for k, v in r.items()}
        if "T" in r:
            r["T"] = int(r["T"])
        out.append(r)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args(argv)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = {k: OUT_DIR / f for k, f in (("runs", "g6-tsw.csv"), ("ev", "g6-tsw-dwells.csv"),
                                         ("sum", "g6-tsw-dwell-summary.csv"),
                                         ("hist", "g6-tsw-rg-hist.csv"), ("man", "g6-tsw.json"))}

    if args.plot_only:
        rows, events, summary = _read(paths["runs"]), _read(paths["ev"]), _read(paths["sum"])
        hrows = list(csv.DictReader(open(paths["hist"])))
        hist = {k: [float(r[k]) for r in hrows] for k in hrows[0] if k != "rg_lo_nm"}
        plot(rows, events, [s for s in summary if s["rule"] == "hyst"], hist,
             json.loads(paths["man"].read_text()), OUT_DIR)
        return 0

    cells = sorted(SWEEP + ENDS)
    with Pool(args.jobs) as pool:
        results = pool.map(both, [name for _, name, _ in cells])
    missing = [name for (_, name, _), r in zip(cells, results) if not r]
    if missing:
        print(f"not finished ({len(missing)}): {', '.join(missing)}")

    rng = np.random.default_rng(SEED)
    rows, events, summary, hist = [], [], [], {}
    for (T, name, kind), r in zip(cells, results):
        if not r:
            continue
        dw = r["dw"]
        for rule, mn in (("hyst", MIN_FRAMES), ("bare", 1)):
            d = np.array([L for L, _ in dw[rule]], dtype=float) * dw["dt_ns"]
            cens = np.array([c for _, c in dw[rule]], dtype=bool)
            events += [{"run": name, "rule": rule, "duration_ns": L, "censored": int(c)}
                       for L, c in zip(d, cens)]
            s = summarise(d, cens, mn, dw["dt_ns"], rng)
            summary.append({"run": name, "T": T, "rule": rule} | s)
            if rule == "hyst":
                hyst = s
        rg = r["rg"]
        for k in ("rgd_sasa_nm2", "rgd_sasa_sem_nm2"):
            rg[k.replace("_nm2", "_per_site_nm2")] = rg.pop(k) / 4
        rows.append({"T": T, "kind": kind} | rg
                    | {k: hyst[k] for k in ("n_dwells", "n_censored", "median_ns",
                                            "km_median_ns", "tau_ns", "tau_lo_ns",
                                            "tau_hi_ns")})
        if kind == "sweep":                 # the heat map is the one-seed sweep only
            hist[str(T)] = r["hist"]

    sw = [x for x in rows if x["kind"] == "sweep"]
    Ts = np.array([x["T"] for x in sw], dtype=float)
    theta = rows[0]["rg_theta_nm"]
    dw0 = next(r for r in results if r)["dw"]
    dt = dw0["dt_ns"]
    manifest = {
        "runs_analysed": [x["run"] for x in rows], "runs_missing": missing,
        "n_residues": dw0["n_residues"], "rg_theta_nm": theta,
        "analysis_window_ns": [DISCARD_NS, DISCARD_NS + dw0["frames"] * dt],
        "frame_ns": dt, "min_dwell_ns": MIN_FRAMES * dt,
        "T_theta_K": crossing(Ts, np.array([x["rg_nm"] for x in sw]), theta),
        "T_half_K": crossing(Ts, np.array([x["collapse_fraction"] for x in sw]), 0.5),
        "primary_rule": {"enter_below": f"{LOW} * Rg_theta", "exit_above": f"{HIGH} * Rg_theta",
                         "min_dwell_frames": MIN_FRAMES},
        "bootstrap": {"resamples": N_BOOT, "seed": SEED},
    }
    _write(paths["runs"], rows)
    _write(paths["ev"], events)
    _write(paths["sum"], summary)
    _write(paths["hist"], [{"rg_lo_nm": f"{BINS[i]:.2f}"} | {t: hist[t][i] for t in hist}
                           for i in range(len(BINS) - 1)])
    paths["man"].write_text(json.dumps(manifest, indent=2))

    print(f"{'run':<16}{'Rg':>6}{'':>6}{'coll':>6}{'n':>5}{'tau':>7}{'tau 95% CI':>14}"
          f"{'RGD':>7}")
    for x in rows:
        print(f"{x['run'][3:]:<16}{x['rg_nm']:>6.2f} ±{x['rg_sem_nm']:.2f}"
              f"{x['collapse_fraction']:>6.0%}{x['n_dwells']:>5}{x['tau_ns']:>7.2f}"
              f"  {x['tau_lo_ns']:>5.2f}-{x['tau_hi_ns']:<5.2f}"
              f"{x['rgd_sasa_per_site_nm2']:>7.2f}")
    print(f"T_theta ≈ {manifest['T_theta_K']:.1f} K, T_half ≈ {manifest['T_half_K']:.1f} K")
    plot(rows, events, [s for s in summary if s["rule"] == "hyst"], hist, manifest, OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
