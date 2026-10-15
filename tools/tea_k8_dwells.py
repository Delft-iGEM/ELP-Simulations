"""Collapsed-state dwells for the 8-lysine arrangements (g6-k8-*) and the original.

    .venv/bin/python tools/tea_k8_dwells.py               # analyse + plot
    .venv/bin/python tools/tea_k8_dwells.py --plot-only   # re-plot from the CSVs

Runs: the ten g6-k8-* arrangements of tools/tea_g6.py (357 aa, 8 VPGKG, only
the lysine positions differ) and g6-orig (the unchanged 4-lysine construct) as
reference, at 285 and 370 K, TEA gamma 6, 400 ns.

Same definitions as tools/tea_lysramp_dwells.py, whose functions this uses:
hysteresis dwells (start below 0.95 Rg_theta, end above 1.05 Rg_theta, minimum
3 frames), a bare-threshold check, right-censoring at the window edges, tau
from censored geometric ML with a 1000-resample bootstrap CI, and the fraction
of frames with Rg < Rg_theta. The x axis is the arrangement - categories, not
a quantity.

Writes, in figures/g6/: g6-k8-dwells.csv, g6-k8-dwell-summary.csv,
g6-k8-dwells.json and g6-k8-dwells-{strip,tau,collapse}.png.
"""

from __future__ import annotations

import argparse
import csv
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from tools.tea_g6 import K8
from tools.tea_lysramp_dwells import (COLD, DISCARD_NS, GRID, HIGH, HOT, INK, LOW,
                                      MIN_FRAMES, MUTED, N_BOOT, SEED, SMALL_N, SURFACE,
                                      TAU_FLOOR, analyse, summarise)

TEMPERATURES = (285, 370)
LABELS = ["orig"] + [f"k8-{k}" for k in K8]
OUT_DIR = Path("figures/g6")


def run_name(label: str, t: int) -> str:
    return f"g6-{label}-{t}k"


def tick(label: str) -> str:
    return "original\n(4 K)" if label == "orig" else label[3:].replace("-", "\n", 1)


# ---- figures ---------------------------------------------------------------------

def _axes(plt, height: float):
    fig, ax = plt.subplots(figsize=(12.0, height))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.grid(color=GRID, linewidth=0.6, axis="y")
    ax.set_axisbelow(True)
    ax.set_xticks(range(len(LABELS)))
    ax.set_xticklabels([tick(lab) for lab in LABELS], fontsize=8.5)
    ax.set_xlim(-0.6, len(LABELS) - 0.4)
    # the reference sits apart from the ten 8-lysine arrangements
    ax.axvline(0.5, color=GRID, linewidth=1.2, zorder=0)
    ax.set_xlabel("Lysine arrangement (8 VPGKG unless noted)", fontsize=9, color=MUTED)
    return fig, ax


def _finish(fig, plt, title: str, caption: str, extra: list, out: Path) -> None:
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, mec=SURFACE, ms=7,
                          label=f"{T} K") for T, c in zip(TEMPERATURES, (COLD, HOT))] + extra
    fig.legend(handles=handles, loc="upper right", ncol=len(handles), frameon=False,
               fontsize=9, labelcolor=INK, bbox_to_anchor=(0.995, 0.995))
    fig.suptitle(title, x=0.01, y=0.995, ha="left", fontsize=11.5, color=INK,
                 fontweight="bold")
    n_lines = caption.count("\n") + 1
    fig.text(0.01, 0.01, caption, fontsize=8, color=MUTED, ha="left", va="bottom",
             linespacing=1.5)
    h = fig.get_size_inches()[1]
    fig.tight_layout(rect=(0, (0.22 * n_lines + 0.1) / h, 1, 1 - 0.45 / h))
    fig.subplots_adjust(top=1 - 0.6 / h)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def plot(events: list[dict], summary: list[dict], m: dict, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = list(zip(TEMPERATURES, (COLD, HOT)))
    off = {TEMPERATURES[0]: -0.17, TEMPERATURES[1]: 0.17}
    by = {(s["label"], s["T"]): s for s in summary}
    rule = (f"Collapsed: Rg < Rg_θ = 0.21·N^0.5 = {m['rg_theta_nm']:.3f} nm "
            f"(N = {m['n_residues']}); dwells with hysteresis (start below {LOW:g}·Rg_θ, "
            f"end above {HIGH:g}·Rg_θ), minimum {MIN_FRAMES} frames "
            f"({m['min_dwell_ns']:.2f} ns); {m['analysis_window_ns'][0]:g}–"
            f"{m['analysis_window_ns'][1]:g} ns analysed, one frame = {m['frame_ns']:.2f} ns.")
    open_marker = plt.Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=MUTED, mew=1.6,
                             ms=7, label=f"n < {SMALL_N} dwells")

    # 1. every dwell, median and IQR
    rng = np.random.default_rng(SEED)
    fig, ax = _axes(plt, 5.6)
    for i, lab in enumerate(LABELS):
        for T, colour in temps:
            s = by.get((lab, T))
            if not s:
                continue
            x = i + off[T]
            d = np.array([e["duration_ns"] for e in events
                          if e["run"] == s["run"] and e["rule"] == "hyst"])
            if len(d):
                ax.scatter(x + rng.uniform(-0.1, 0.1, len(d)), d, s=6, color=colour,
                           alpha=0.25, linewidths=0, zorder=2)
            if s["n_dwells"]:
                small = s["n_dwells"] < SMALL_N
                ax.plot([x, x], [s["q25_ns"], s["q75_ns"]], color=colour, lw=2.2, zorder=3,
                        solid_capstyle="butt")
                ax.plot(x, s["median_ns"], "o", ms=7.5, zorder=4,
                        mfc=SURFACE if small else colour, mec=colour, mew=1.8)
    ax.set_yscale("log")
    ax.set_ylabel("Collapsed dwell (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, "Collapsed-state dwells, 8-lysine arrangements",
            rule + "\nOne faint point per dwell, jittered; solid marker = median, "
                   "bar = interquartile range.",
            [plt.Line2D([], [], marker="o", ls="", mfc=MUTED, mec=MUTED, ms=7,
                        label="median, bar = IQR"), open_marker],
            out_dir / "g6-k8-dwells-strip.png")

    # 2. tau, each with its CI bar
    fig, ax = _axes(plt, 4.8)
    dt = m["frame_ns"]
    for i, lab in enumerate(LABELS):
        for T, colour in temps:
            s = by.get((lab, T))
            if not s or not np.isfinite(s["tau_ns"]):
                continue
            x = i + off[T]
            small = s["n_dwells"] < SMALL_N
            ax.plot([x, x], [max(s["tau_lo_ns"], TAU_FLOOR), max(s["tau_hi_ns"], TAU_FLOOR)],
                    color=colour, lw=1.4, zorder=2, solid_capstyle="butt")
            if s["tau_ns"] < TAU_FLOOR:
                ax.plot(x, TAU_FLOOR, "v", ms=6.5, zorder=3, clip_on=False,
                        mfc=SURFACE if small else colour, mec=colour, mew=1.4)
            else:
                ax.plot(x, s["tau_ns"], "o", ms=6.5, zorder=3,
                        mfc=SURFACE if small else colour, mec=colour, mew=1.6)
    ax.set_yscale("log")
    ax.set_ylim(TAU_FLOOR, 3.0)
    ax.axhspan(TAU_FLOOR, dt, color=GRID, alpha=0.6, linewidth=0, zorder=0)
    ax.text(0.01, TAU_FLOOR * 1.25, f" below one frame ({dt:.2f} ns): not resolved",
            transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8,
            color=MUTED)
    ax.set_ylabel("Dwell time τ (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, "Collapsed-state dwell time τ, 8-lysine arrangements",
            rule + "\nτ: exponential (geometric ML on whole frames above the minimum, dwells "
                   "touching the window edges right-censored); bar = 95 % bootstrap CI "
                   f"({N_BOOT} resamples);\n▼ = τ of 0 (every dwell at the minimum).",
            [open_marker], out_dir / "g6-k8-dwells-tau.png")

    # 3. fraction of time collapsed
    fig, ax = _axes(plt, 4.4)
    for T, colour in temps:
        pts = [(i + off[T], by[(lab, T)]["collapse_fraction"])
               for i, lab in enumerate(LABELS) if (lab, T) in by]
        ax.plot(*zip(*pts), "o", ms=7.5, color=colour, mec=SURFACE, mew=1.3, zorder=3)
    ax.set_ylim(-0.03, 1.0)
    ax.set_ylabel("Fraction of time collapsed", fontsize=9.5, color=INK)
    _finish(fig, plt, "Fraction of time collapsed, 8-lysine arrangements",
            f"Fraction of the {m['analysis_window_ns'][0]:g}–{m['analysis_window_ns'][1]:g} ns "
            f"frames with Rg < Rg_θ = {m['rg_theta_nm']:.3f} nm (no hysteresis, no minimum "
            "length).", [], out_dir / "g6-k8-dwells-collapse.png")


# ---- main ------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args(argv)
    ev_csv = OUT_DIR / "g6-k8-dwells.csv"
    sum_csv = OUT_DIR / "g6-k8-dwell-summary.csv"
    man_json = OUT_DIR / "g6-k8-dwells.json"

    if args.plot_only:
        events = [dict(r, duration_ns=float(r["duration_ns"])) for r in csv.DictReader(open(ev_csv))]
        summary = [{k: (v if k in ("run", "label", "rule") else
                        int(v) if k in ("T", "n_dwells", "n_censored") else float(v))
                    for k, v in r.items()} for r in csv.DictReader(open(sum_csv))]
        plot(events, [s for s in summary if s["rule"] == "hyst"],
             json.loads(man_json.read_text()), OUT_DIR)
        return 0

    cells = [(lab, T) for lab in LABELS for T in TEMPERATURES]
    with Pool(args.jobs) as pool:
        results = pool.map(analyse, [run_name(lab, T) for lab, T in cells])
    missing = [run_name(*c) for c, r in zip(cells, results) if not r]
    if missing:
        print(f"not finished ({len(missing)}): {', '.join(missing)}")
    done = [(c, r) for c, r in zip(cells, results) if r]

    rng = np.random.default_rng(SEED)
    events, summary = [], []
    for (lab, T), r in done:
        for rule, mn in (("hyst", MIN_FRAMES), ("bare", 1)):
            d = np.array([L for L, _ in r[rule]], dtype=float) * r["dt_ns"]
            cens = np.array([c for _, c in r[rule]], dtype=bool)
            events += [{"run": r["run"], "rule": rule, "duration_ns": L, "censored": int(c)}
                       for L, c in zip(d, cens)]
            summary.append({"run": r["run"], "label": lab, "T": T, "rule": rule,
                            "collapse_fraction": r["collapse_fraction"]}
                           | summarise(d, cens, mn, r["dt_ns"], rng))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for path, rows in ((ev_csv, events), (sum_csv, summary)):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    r0 = done[0][1]
    dt = r0["dt_ns"]
    manifest = {
        "runs_analysed": [r["run"] for _, r in done], "runs_missing": missing,
        "n_residues": r0["n_residues"], "rg_theta_nm": r0["theta_nm"],
        "rg_theta_rule": "0.21 * N**0.5 nm",
        "analysis_window_ns": [DISCARD_NS, DISCARD_NS + r0["frames"] * dt],
        "frame_ns": dt,
        "primary_rule": {"name": "hyst", "enter_below": f"{LOW} * Rg_theta",
                         "exit_above": f"{HIGH} * Rg_theta", "min_dwell_frames": MIN_FRAMES},
        "check_rule": {"name": "bare", "collapsed_if": "Rg < Rg_theta", "min_dwell_frames": 1},
        "min_dwell_ns": MIN_FRAMES * dt,
        "censoring": "dwells touching the first or last analysed frame are right-censored; "
                     "tau from censored geometric ML on frames above the minimum, "
                     "p = sum(L - L_min) / (sum(L - L_min) + n_completed), tau = -dt / ln p; "
                     "Kaplan-Meier median also reported",
        "bootstrap": {"resamples": N_BOOT, "ci": "percentile 2.5-97.5", "seed": SEED},
        "small_n_flag": SMALL_N,
    }
    man_json.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {ev_csv}, {sum_csv}, {man_json}")

    print(f"{'run':<26}{'rule':>5}{'n':>5}{'cens':>5}{'median':>8}{'IQR':>14}"
          f"{'KM med':>8}{'tau':>7}{'tau 95% CI':>15}{'coll':>6}")
    for s in summary:
        print(f"{s['run'][3:]:<26}{s['rule']:>5}{s['n_dwells']:>5}{s['n_censored']:>5}"
              f"{s['median_ns']:>8.2f}  {s['q25_ns']:>5.2f}-{s['q75_ns']:<6.2f}"
              f"{s['km_median_ns']:>8.2f}{s['tau_ns']:>7.2f}"
              f"  {s['tau_lo_ns']:>5.2f}-{s['tau_hi_ns']:<6.2f}{s['collapse_fraction']:>6.0%}")
    plot(events, [s for s in summary if s["rule"] == "hyst"], manifest, OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
