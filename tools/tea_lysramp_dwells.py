"""Collapsed-state dwell times and collapse fraction for the lysine ramp.

    .venv/bin/python tools/tea_lysramp_dwells.py               # analyse + plot
    .venv/bin/python tools/tea_lysramp_dwells.py --plot-only   # re-plot from the CSVs

Runs: tools/tea_lysramp.py (g6-kspc*/g6-kend*, 357 aa, 290 / 370 K, 400 ns).
Frames before 50 ns are discarded, as in every g6 analysis.

Collapsed state
---------------
Rg_theta = 0.21 N^0.5 nm (tools/tea_analysis.rg_theta), N = 357 -> 3.968 nm.

  hysteresis (primary)  a dwell starts when Rg drops below 0.95 Rg_theta and
                        ends when Rg rises above 1.05 Rg_theta; dwells shorter
                        than MIN_FRAMES frames are discarded.
  bare (check)          a dwell is any run of consecutive frames with
                        Rg < Rg_theta, no minimum length.

A dwell that is in progress at the first analysed frame or the last frame is
right-censored: its true length is longer than observed.

Summaries per (series, n_lys, T)
--------------------------------
  median, IQR       of the observed dwell lengths, censored ones included as
                    observed (a lower bound - they are few, see n_censored)
  km_median         Kaplan-Meier median, censored dwells as censored
  tau               exponential dwell time. Dwells are whole frames, so the
                    excess over the minimum, L - L_min, is fitted as geometric
                    (the discrete exponential; memoryless, so the cut at L_min
                    does not bias it), maximum likelihood with right censoring:
                        p = sum(L - L_min) / (sum(L - L_min) + n_completed)
                        tau = -dt / ln p
                    95 % CI: percentile bootstrap, 1000 resamples of the
                    group's dwells. tau below one frame (dt) is not resolved.
  collapse_fraction fraction of analysed frames with Rg < Rg_theta (bare rule,
                    the same number as g6-lysramp.csv)

Writes, in figures/g6-kspc-kend/: g6-lysramp-dwells.csv (one row per dwell),
g6-lysramp-dwell-summary.csv, g6-lysramp-dwells.json (manifest), and the figures
g6-lysramp-dwells-{strip,tau,collapse}.png.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_analysis import rg_theta
from tools.tea_lysramp import COUNTS, FREE, SERIES, TEMPERATURES

DISCARD_NS = 50.0
LOW, HIGH = 0.95, 1.05
MIN_FRAMES = 3
N_BOOT = 1000
SMALL_N = 10
SEED = 20261014


def run_name(series: str, n: int, t: float) -> str:
    return f"g6-{series}{n:02d}-{int(t)}k"


def rg_series(name: str) -> tuple[np.ndarray, float, int] | None:
    """Per-frame Rg (nm) after DISCARD_NS, the frame time (ns), and N."""
    import mdtraj as md

    from tools.analysis_batch import trajectory_paths
    from tools.metadata import read_metadata
    from tools.tea_g6_analysis import FRAMES_EXPECTED, _whole

    try:
        top, dcd = trajectory_paths(name)
    except FileNotFoundError:
        return None
    if not Path(dcd).is_file():
        return None
    traj = md.load(str(dcd), top=str(top))
    if traj.n_frames < FRAMES_EXPECTED:
        return None
    meta = read_metadata(Path("simulations") / name / "runtime")
    dt = float(meta["time_per_frame_ps"]) / 1000.0
    traj = traj[np.arange(traj.n_frames) * dt >= DISCARD_NS]
    box = traj.unitcell_lengths[0]
    xyz = np.stack([_whole(f, box) for f in traj.xyz])
    dev = xyz - xyz.mean(axis=1, keepdims=True)
    return np.sqrt((dev ** 2).sum(axis=2).mean(axis=1)), dt, traj.n_atoms


def collapsed_mask(rg: np.ndarray, theta: float, hysteresis: bool) -> np.ndarray:
    if not hysteresis:
        return rg < theta
    mask = np.empty(len(rg), dtype=bool)
    # In the band at the first frame, take the side of theta it is on.
    state = bool(rg[0] < theta)
    for i, r in enumerate(rg):
        if state and r > HIGH * theta:
            state = False
        elif not state and r < LOW * theta:
            state = True
        mask[i] = state
    return mask


def dwells(mask: np.ndarray, min_frames: int) -> list[tuple[int, bool]]:
    """(length in frames, censored) for each run of True."""
    edges = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
    starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    return [(int(e - s), bool(s == 0 or e == len(mask)))
            for s, e in zip(starts, ends) if e - s >= min_frames]


def analyse(name: str) -> dict | None:
    got = rg_series(name)
    if got is None:
        return None
    rg, dt, n = got
    theta = rg_theta(n)
    out = {"run": name, "dt_ns": dt, "n_residues": n, "theta_nm": theta,
           "frames": len(rg), "collapse_fraction": float((rg < theta).mean())}
    for rule, hyst, mn in (("hyst", True, MIN_FRAMES), ("bare", False, 1)):
        out[rule] = dwells(collapsed_mask(rg, theta, hyst), mn)
    return out


# ---- statistics ----------------------------------------------------------------

def tau_mle(frames: np.ndarray, cens: np.ndarray, min_frames: int, dt: float) -> float:
    """Censored geometric ML for the excess frames, as a time constant in ns."""
    k = int((~cens).sum())
    if not k:
        return float("nan")
    excess = float((frames - min_frames).sum())
    p = excess / (excess + k)
    return -dt / math.log(p) if p > 0 else 0.0


def km_median(d: np.ndarray, cens: np.ndarray) -> float:
    s = 1.0
    for t in np.unique(d):
        at = d >= t
        events = int(((d == t) & ~cens).sum())
        if events:
            s *= 1 - events / int(at.sum())
            if s <= 0.5:
                return float(t)
    return float("nan")                       # survival never reaches 0.5


def summarise(d: np.ndarray, cens: np.ndarray, min_frames: int, dt: float,
              rng: np.random.Generator) -> dict:
    row = {"n_dwells": len(d), "n_censored": int(cens.sum())}
    if not len(d):
        return row | {k: float("nan") for k in
                      ("median_ns", "q25_ns", "q75_ns", "km_median_ns",
                       "tau_ns", "tau_lo_ns", "tau_hi_ns")}
    q25, med, q75 = np.percentile(d, [25, 50, 75])
    frames = np.rint(d / dt).astype(int)
    boot = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(d), len(d))
        boot.append(tau_mle(frames[i], cens[i], min_frames, dt))
    boot = np.array(boot)
    boot = boot[np.isfinite(boot)]
    lo, hi = (np.percentile(boot, [2.5, 97.5]) if len(boot) > N_BOOT // 2
              else (float("nan"), float("nan")))
    return row | {"median_ns": float(med), "q25_ns": float(q25), "q75_ns": float(q75),
                  "km_median_ns": km_median(d, cens),
                  "tau_ns": tau_mle(frames, cens, min_frames, dt),
                  "tau_lo_ns": float(lo), "tau_hi_ns": float(hi)}


def group_key(name: str) -> tuple[str, int, int]:
    return name[3:7], int(name[7:9]), int(name[10:13])


# ---- figure ----------------------------------------------------------------------

COLD = "#4269d0"
HOT = "#c1442e"
INK = "#1a1a1a"
MUTED = "#6b6f76"
GRID = "#dcdee2"
SURFACE = "#fcfcfb"
TAU_FLOOR = 0.005
TITLE = {"kspc": "Lysines spaced evenly", "kend": "Lysines filled in from both ends"}


def _axes(plt, n_rows_ratio: float = 1.0):
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6 * n_rows_ratio), sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, series in zip(axes, SERIES):
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(TITLE[series], fontsize=10.5, color=INK, fontweight="bold")
        ax.set_xlabel("Lysine fraction (VPGKG / 49 replaceable pentads)",
                      fontsize=9, color=MUTED)
        ax.set_xticks([n / len(FREE) for n in COUNTS])
        ax.set_xticklabels([f"{n / len(FREE):.2f}" for n in COUNTS])
        ax.set_xlim(0.02, 1.06)
    return fig, axes


def _finish(fig, plt, title: str, caption: str, extra: list, out: Path) -> None:
    temps = [(int(T), c) for T, c in zip(TEMPERATURES, (COLD, HOT))]
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, mec=SURFACE, ms=7,
                          label=f"{T} K") for T, c in temps] + extra
    fig.legend(handles=handles, loc="upper right", ncol=len(handles), frameon=False,
               fontsize=9, labelcolor=INK, bbox_to_anchor=(0.995, 0.995))
    fig.suptitle(title, x=0.01, y=0.995, ha="left", fontsize=11.5, color=INK,
                 fontweight="bold")
    n_lines = caption.count("\n") + 1
    fig.text(0.01, 0.01, caption, fontsize=8, color=MUTED, ha="left", va="bottom",
             linespacing=1.5)
    h = fig.get_size_inches()[1]
    fig.tight_layout(rect=(0, (0.22 * n_lines + 0.1) / h, 1, 1 - 0.45 / h))
    top = 1 - 0.75 / h                       # tight_layout leaves room for the legend twice
    fig.subplots_adjust(top=top)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def plot(events: list[dict], summary: list[dict], manifest: dict, out_dir: Path,
         prefix: str) -> None:
    """Three figures: <prefix>s-strip.png, <prefix>s-tau.png, <prefix>s-collapse.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = [(int(T), c) for T, c in zip(TEMPERATURES, (COLD, HOT))]
    off = {temps[0][0]: -0.018, temps[1][0]: 0.018}
    m = manifest
    rule = (f"Collapsed: Rg < Rg_θ = 0.21·N^0.5 = {m['rg_theta_nm']:.3f} nm "
            f"(N = {m['n_residues']}); dwells with hysteresis (start below {LOW:g}·Rg_θ, "
            f"end above {HIGH:g}·Rg_θ), minimum {MIN_FRAMES} frames "
            f"({m['min_dwell_ns']:.2f} ns); {m['analysis_window_ns'][0]:g}–"
            f"{m['analysis_window_ns'][1]:g} ns analysed, one frame = {m['frame_ns']:.2f} ns.")
    open_marker = plt.Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=MUTED, mew=1.6,
                             ms=7, label=f"n < {SMALL_N} dwells")
    groups = {}
    for s in summary:
        groups.setdefault((s["series"], s["T"]), []).append(s)
    for v in groups.values():
        v.sort(key=lambda s: s["n_lys"])

    # 1. every dwell, median and IQR
    rng = np.random.default_rng(SEED)
    fig, axes = _axes(plt, 1.25)
    for ax, series in zip(axes, SERIES):
        for T, colour in temps:
            for s in groups.get((series, T), []):
                x = s["lys_fraction"] + off[T]
                d = np.array([e["duration_ns"] for e in events
                              if e["run"] == s["run"] and e["rule"] == "hyst"])
                if len(d):
                    ax.scatter(x + rng.uniform(-0.011, 0.011, len(d)), d, s=6,
                               color=colour, alpha=0.25, linewidths=0, zorder=2)
                if s["n_dwells"]:
                    small = s["n_dwells"] < SMALL_N
                    ax.plot([x, x], [s["q25_ns"], s["q75_ns"]], color=colour, lw=2.2,
                            zorder=3, solid_capstyle="butt")
                    ax.plot(x, s["median_ns"], "o", ms=7.5, zorder=4,
                            mfc=SURFACE if small else colour, mec=colour, mew=1.8)
        ax.set_yscale("log")
    axes[0].set_ylabel("Collapsed dwell (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, "Collapsed-state dwells, lysine ramp on the 357-aa RGD construct",
            rule + "\nOne faint point per dwell, jittered; solid marker = median, "
                   "bar = interquartile range.",
            [plt.Line2D([], [], marker="o", ls="", mfc=MUTED, mec=MUTED, ms=7,
                        label="median, bar = IQR"), open_marker],
            out_dir / f"{prefix}s-strip.png")

    # 2. tau, each point with its own CI bar, no connecting lines
    fig, axes = _axes(plt)
    dt = m["frame_ns"]
    for ax, series in zip(axes, SERIES):
        for T, colour in temps:
            for s in groups.get((series, T), []):
                if not np.isfinite(s["tau_ns"]):
                    continue
                x = s["lys_fraction"] + off[T]
                small = s["n_dwells"] < SMALL_N
                ax.plot([x, x], [max(s["tau_lo_ns"], TAU_FLOOR),
                                 max(s["tau_hi_ns"], TAU_FLOOR)],
                        color=colour, lw=1.4, zorder=2, solid_capstyle="butt")
                if s["tau_ns"] < TAU_FLOOR:          # zero: pinned to the axis floor
                    ax.plot(x, TAU_FLOOR, "v", ms=6.5, zorder=3, clip_on=False,
                            mfc=SURFACE if small else colour, mec=colour, mew=1.4)
                else:
                    ax.plot(x, s["tau_ns"], "o", ms=6.5, zorder=3,
                            mfc=SURFACE if small else colour, mec=colour, mew=1.6)
        ax.set_yscale("log")
        ax.set_ylim(TAU_FLOOR, 2.0)
        ax.axhspan(TAU_FLOOR, dt, color=GRID, alpha=0.6, linewidth=0, zorder=0)
        ax.text(0.01, TAU_FLOOR * 1.25, f" below one frame ({dt:.2f} ns): not resolved",
                transform=ax.get_yaxis_transform(), ha="left", va="bottom",
                fontsize=8, color=MUTED)
    axes[0].set_ylabel("Dwell time τ (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, "Collapsed-state dwell time τ, lysine ramp on the 357-aa RGD construct",
            rule + "\nτ: exponential (geometric ML on whole frames above the minimum, dwells "
                   "touching the window edges right-censored); bar = 95 % bootstrap CI "
                   f"({N_BOOT} resamples);\n▼ = τ of 0 (every dwell at the minimum).",
            [open_marker], out_dir / f"{prefix}s-tau.png")

    # 3. fraction of time collapsed
    fig, axes = _axes(plt, 0.9)
    for ax, series in zip(axes, SERIES):
        for T, colour in temps:
            rows = groups.get((series, T), [])
            ax.plot([s["lys_fraction"] + off[T] for s in rows],
                    [s["collapse_fraction"] for s in rows], "o", ms=7.5,
                    color=colour, mec=SURFACE, mew=1.3, zorder=3)
        ax.set_ylim(-0.03, 1.0)
    axes[0].set_ylabel("Fraction of time collapsed", fontsize=9.5, color=INK)
    _finish(fig, plt, "Fraction of time collapsed, lysine ramp on the 357-aa RGD construct",
            f"Fraction of the {m['analysis_window_ns'][0]:g}–{m['analysis_window_ns'][1]:g} ns "
            f"frames with Rg < Rg_θ = {m['rg_theta_nm']:.3f} nm (no hysteresis, no minimum "
            "length).", [], out_dir / f"{prefix}s-collapse.png")


# ---- main -------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prefix", default="g6-lysramp-dwell")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args(argv)
    out_dir = Path("figures/g6-kspc-kend")
    ev_csv = out_dir / f"{args.prefix}s.csv"
    sum_csv = out_dir / f"{args.prefix}-summary.csv"
    man_json = out_dir / f"{args.prefix}s.json"

    if args.plot_only:
        events = [dict(r, duration_ns=float(r["duration_ns"])) for r in csv.DictReader(open(ev_csv))]
        summary = []
        for r in csv.DictReader(open(sum_csv)):
            summary.append({k: (v if k in ("run", "series", "rule") else
                                int(v) if k in ("n_lys", "T", "n_dwells", "n_censored")
                                else float(v)) for k, v in r.items()})
        plot(events, [s for s in summary if s["rule"] == "hyst"],
             json.loads(man_json.read_text()), out_dir, args.prefix)
        return 0

    names = [run_name(s, n, t) for s in SERIES for n in COUNTS for t in TEMPERATURES]
    with Pool(args.jobs) as pool:
        results = pool.map(analyse, names)
    done = [r for r in results if r]
    missing = [n for n, r in zip(names, results) if not r]
    if missing:
        print(f"not finished ({len(missing)}): {', '.join(missing)}")

    rng = np.random.default_rng(SEED)
    events, summary = [], []
    for r in done:
        series, n_lys, T = group_key(r["run"])
        for rule, mn in (("hyst", MIN_FRAMES), ("bare", 1)):
            d = np.array([L for L, _ in r[rule]], dtype=float) * r["dt_ns"]
            cens = np.array([c for _, c in r[rule]], dtype=bool)
            for L, c in zip(d, cens):
                events.append({"run": r["run"], "rule": rule, "duration_ns": L,
                               "censored": int(c)})
            summary.append({"run": r["run"], "series": series, "n_lys": n_lys, "T": T,
                            "lys_fraction": n_lys / len(FREE), "rule": rule,
                            "collapse_fraction": r["collapse_fraction"]}
                           | summarise(d, cens, mn, r["dt_ns"], rng))

    with open(ev_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(events[0]))
        w.writeheader()
        w.writerows(events)
    with open(sum_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    dt = done[0]["dt_ns"]
    manifest = {
        "runs_analysed": [r["run"] for r in done], "runs_missing": missing,
        "n_residues": done[0]["n_residues"], "rg_theta_nm": done[0]["theta_nm"],
        "rg_theta_rule": "0.21 * N**0.5 nm",
        "analysis_window_ns": [DISCARD_NS, DISCARD_NS + done[0]["frames"] * dt],
        "frame_ns": dt,
        "primary_rule": {"name": "hyst", "enter_below": f"{LOW} * Rg_theta",
                         "exit_above": f"{HIGH} * Rg_theta",
                         "min_dwell_frames": MIN_FRAMES},
        "check_rule": {"name": "bare", "collapsed_if": "Rg < Rg_theta",
                       "min_dwell_frames": 1},
        "min_dwell_ns": MIN_FRAMES * dt,
        "censoring": "dwells touching the first or last analysed frame are right-"
                     "censored; tau from censored geometric ML on frames above the "
                     "minimum, p = sum(L - L_min) / (sum(L - L_min) + n_completed), "
                     "tau = -dt / ln p; Kaplan-Meier median also reported",
        "bootstrap": {"resamples": N_BOOT, "ci": "percentile 2.5-97.5", "seed": SEED},
        "small_n_flag": SMALL_N,
    }
    man_json.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {ev_csv}, {sum_csv}, {man_json}")

    hdr = (f"{'run':<16}{'rule':>5}{'n':>5}{'cens':>5}{'median':>8}{'IQR':>14}"
           f"{'KM med':>8}{'tau':>7}{'tau 95% CI':>15}{'coll':>6}")
    print(hdr)
    for s in sorted(summary, key=lambda s: (s["series"], s["T"], s["n_lys"], s["rule"] != "hyst")):
        print(f"{s['run'][3:]:<16}{s['rule']:>5}{s['n_dwells']:>5}{s['n_censored']:>5}"
              f"{s['median_ns']:>8.2f}  {s['q25_ns']:>5.2f}-{s['q75_ns']:<6.2f}"
              f"{s['km_median_ns']:>8.2f}{s['tau_ns']:>7.2f}"
              f"  {s['tau_lo_ns']:>5.2f}-{s['tau_hi_ns']:<6.2f}{s['collapse_fraction']:>6.0%}")
    plot(events, [s for s in summary if s["rule"] == "hyst"], manifest, out_dir, args.prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
