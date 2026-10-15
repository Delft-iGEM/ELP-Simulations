"""Rg and collapsed-state dwells for the random arrangements (tools/tea_random.py).

    .venv/bin/python tools/tea_random_analysis.py               # analyse + plot
    .venv/bin/python tools/tea_random_analysis.py --plot-only   # re-plot from the CSVs

Runs, all 357 aa, 285 and 370 K, TEA gamma 6, 400 ns, 50 ns discarded:

  orig    g6-orig, the unchanged construct (4 K)                     1 sequence
  k8      g6-k8-*, the ten designed 8-lysine arrangements            10
  k8r     g6-k8r*, 8 lysines at random pentads, RGD fixed            10
  rk      g6-rk*, the 4 RGD blocks and 4 lysines both placed at random 10

The designed and random k8 sets are re-analysed here so every number in the
figures comes from the same code. Rg and its block SEM come from
tools/tea_g6_analysis.analyse_run; dwells, tau and collapse fraction from
tools/tea_lysramp_dwells (hysteresis 0.95 / 1.05 Rg_theta, minimum 3 frames,
censored geometric tau with a 1000-resample bootstrap CI).

Writes, in figures/g6-random/: g6-random.csv (one row per run: Rg, SASA, tau,
collapse fraction), g6-random-dwells.csv (one row per dwell),
g6-random-dwell-summary.csv, g6-random-dwells.json (manifest) and the figures
g6-random-{rg,collapse,tau,dwells-strip}.png.
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
GROUPS = {
    "orig": ("original\n(4 K)", ["orig"]),
    "k8": ("8 K, designed\npositions", [f"k8-{k}" for k in K8]),
    "k8r": ("8 K, random\npositions", [f"k8r{i:02d}" for i in range(1, 11)]),
    "rk": ("4 K + 4 RGD,\nrandom positions", [f"rk{i:02d}" for i in range(1, 11)]),
}
OUT_DIR = Path("figures/g6-random")
T_OFF = {TEMPERATURES[0]: -0.2, TEMPERATURES[1]: 0.2}
SPREAD = 0.13                     # arrangements fan out over +-SPREAD inside a group


def run_name(label: str, t: int) -> str:
    return f"g6-{label}-{t}k"


def both(name: str) -> dict | None:
    from tools.tea_g6_analysis import analyse_run

    row = analyse_run(name)
    if row is None:
        return None
    return {"rg": row, "dw": analyse(name)}


# ---- figures ---------------------------------------------------------------------

def _style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.grid(color=GRID, linewidth=0.6, axis="y")
    ax.set_axisbelow(True)


def _group_axes(plt, height: float):
    fig, ax = plt.subplots(figsize=(8.4, height))
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels([g[0] for g in GROUPS.values()], fontsize=8.5)
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    for x in np.arange(0.5, len(GROUPS) - 1):
        ax.axvline(x, color=GRID, linewidth=1.2, zorder=0)
    return fig, ax


def _x(group: str, label: str, T: int) -> float:
    labels = GROUPS[group][1]
    g = list(GROUPS).index(group)
    i = labels.index(label)
    fan = 0.0 if len(labels) == 1 else -SPREAD + 2 * SPREAD * i / (len(labels) - 1)
    return g + T_OFF[T] + fan


def _finish(fig, plt, title: str, caption: str, extra: list, out: Path) -> None:
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, mec=SURFACE, ms=7,
                          label=f"{T} K") for T, c in zip(TEMPERATURES, (COLD, HOT))] + extra
    fig.legend(handles=handles, loc="upper right", ncol=len(handles), frameon=False,
               fontsize=8.5, labelcolor=INK, bbox_to_anchor=(0.995, 0.93))
    fig.suptitle(title, x=0.01, y=0.995, ha="left", fontsize=11.5, color=INK,
                 fontweight="bold")
    n_lines = caption.count("\n") + 1
    fig.text(0.01, 0.01, caption, fontsize=7.8, color=MUTED, ha="left", va="bottom",
             linespacing=1.5)
    h = fig.get_size_inches()[1]
    fig.tight_layout(rect=(0, (0.21 * n_lines + 0.1) / h, 1, 1 - 0.45 / h))
    fig.subplots_adjust(top=1 - 0.75 / h)
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def _median_ticks(ax, rows: list[dict], key: str):
    """A short bar at each group's median, where at least 5 sequences have a value."""
    for group, (_, labels) in GROUPS.items():
        if len(labels) == 1:
            continue
        for T in TEMPERATURES:
            v = [r[key] for r in rows if r["group"] == group and r["T"] == T
                 and np.isfinite(r[key])]
            if len(v) >= 5:
                c = list(GROUPS).index(group) + T_OFF[T]
                ax.plot([c - SPREAD - 0.04, c + SPREAD + 0.04], [np.median(v)] * 2,
                        color=INK, lw=1.3, zorder=4, solid_capstyle="butt")


def plot(rows: list[dict], events: list[dict], summary: list[dict], m: dict,
         out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temps = list(zip(TEMPERATURES, (COLD, HOT)))
    colour = dict(temps)
    out_dir.mkdir(parents=True, exist_ok=True)
    med = plt.Line2D([], [], color=INK, lw=1.3, label="group median")
    open_marker = plt.Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=MUTED, mew=1.6,
                             ms=7, label=f"n < {SMALL_N} dwells")
    window = f"{m['analysis_window_ns'][0]:g}–{m['analysis_window_ns'][1]:g} ns"
    rule = (f"Collapsed: Rg < Rg_θ = 0.21·N^0.5 = {m['rg_theta_nm']:.3f} nm "
            f"(N = {m['n_residues']});\ndwells with hysteresis (start below {LOW:g}·Rg_θ, "
            f"end above {HIGH:g}·Rg_θ), minimum {MIN_FRAMES} frames "
            f"({m['min_dwell_ns']:.2f} ns); {window} analysed.")
    one = "One point per sequence; the original sits alone at the left."

    # 1. Rg
    fig, ax = _group_axes(plt, 4.8)
    for r in rows:
        x = _x(r["group"], r["label"], r["T"])
        ax.errorbar(x, r["rg_nm"], yerr=r["rg_sem_nm"], fmt="o", ms=6.5,
                    color=colour[r["T"]], elinewidth=1, capsize=0, mec=SURFACE, mew=1.2,
                    zorder=3)
    _median_ticks(ax, rows, "rg_nm")
    theta = m["rg_theta_nm"]
    ax.axhline(theta, color=MUTED, linewidth=1, linestyle=(0, (1, 2)), zorder=1)
    ax.text(0.995, theta, f"θ-state Rg {theta:.2f} nm", transform=ax.get_yaxis_transform(),
            va="bottom", ha="right", fontsize=8, color=MUTED)
    ax.set_ylabel("Radius of gyration (nm)", fontsize=9.5, color=INK)
    _finish(fig, plt, "Radius of gyration, designed vs random arrangements",
            f"Mean Rg over {window}; bar = SEM from 20 block averages. " + one,
            [med], out_dir / "g6-random-rg.png")

    # 2. collapse fraction
    fig, ax = _group_axes(plt, 4.4)
    for r in rows:
        ax.plot(_x(r["group"], r["label"], r["T"]), r["collapse_fraction"], "o", ms=6.5,
                color=colour[r["T"]], mec=SURFACE, mew=1.2, zorder=3)
    _median_ticks(ax, rows, "collapse_fraction")
    ax.set_ylim(-0.03, 1.0)
    ax.set_ylabel("Fraction of time collapsed", fontsize=9.5, color=INK)
    _finish(fig, plt, "Fraction of time collapsed, designed vs random arrangements",
            f"Fraction of the {window} frames with Rg < Rg_θ = {theta:.3f} nm (no hysteresis, "
            "no minimum length).\n" + one, [med], out_dir / "g6-random-collapse.png")

    # 3. tau
    fig, ax = _group_axes(plt, 4.8)
    dt = m["frame_ns"]
    for r in rows:
        if not np.isfinite(r["tau_ns"]):
            continue
        x, c = _x(r["group"], r["label"], r["T"]), colour[r["T"]]
        small = r["n_dwells"] < SMALL_N
        ax.plot([x, x], [max(r["tau_lo_ns"], TAU_FLOOR), max(r["tau_hi_ns"], TAU_FLOOR)],
                color=c, lw=1.2, zorder=2, solid_capstyle="butt")
        if r["tau_ns"] < TAU_FLOOR:
            ax.plot(x, TAU_FLOOR, "v", ms=6, zorder=3, clip_on=False,
                    mfc=SURFACE if small else c, mec=c, mew=1.4)
        else:
            ax.plot(x, r["tau_ns"], "o", ms=6, zorder=3,
                    mfc=SURFACE if small else c, mec=c, mew=1.5)
    _median_ticks(ax, [r for r in rows if r["n_dwells"] >= SMALL_N], "tau_ns")
    ax.set_yscale("log")
    ax.set_ylim(TAU_FLOOR, 3.0)
    ax.axhspan(TAU_FLOOR, dt, color=GRID, alpha=0.6, linewidth=0, zorder=0)
    ax.text(0.01, TAU_FLOOR * 1.25, f" below one frame ({dt:.2f} ns): not resolved",
            transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=8,
            color=MUTED)
    ax.set_ylabel("Dwell time τ (ns)", fontsize=9.5, color=INK)
    _finish(fig, plt, "Collapsed-state dwell time τ, designed vs random arrangements",
            rule + "\nτ: exponential (censored geometric ML); bar = 95 % bootstrap CI "
                   f"({N_BOOT} resamples); ▼ = τ of 0. Medians only over sequences with "
                   f"≥ {SMALL_N} dwells\n(drawn only where ≥ 5 sequences qualify). " + one,
            [med, open_marker], out_dir / "g6-random-tau.png")

    # 4. every dwell of every run, per sequence
    labels = [lab for g, (_, labs) in GROUPS.items() for lab in labs]
    by = {(s["label"], s["T"]): s for s in summary}
    rng = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(15.0, 5.6))
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([lab if lab == "orig" else lab.replace("k8-", "").split("-")[0]
                        for lab in labels], fontsize=7.5, rotation=90)
    ax.set_xlim(-0.6, len(labels) - 0.4)
    start = 0
    for g, (name, labs) in GROUPS.items():
        if start:
            ax.axvline(start - 0.5, color=GRID, linewidth=1.2, zorder=0)
        ax.text(start + (len(labs) - 1) / 2, 1.0, name.replace("\n", " "),
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8.5,
                color=MUTED)
        start += len(labs)
    for i, lab in enumerate(labels):
        for T, c in temps:
            s = by.get((lab, T))
            if not s:
                continue
            x = i + T_OFF[T] * 0.9
            d = np.array([e["duration_ns"] for e in events
                          if e["run"] == s["run"] and e["rule"] == "hyst"])
            if len(d):
                ax.scatter(x + rng.uniform(-0.08, 0.08, len(d)), d, s=5, color=c,
                           alpha=0.25, linewidths=0, zorder=2)
            if s["n_dwells"]:
                small = s["n_dwells"] < SMALL_N
                ax.plot([x, x], [s["q25_ns"], s["q75_ns"]], color=c, lw=2, zorder=3,
                        solid_capstyle="butt")
                ax.plot(x, s["median_ns"], "o", ms=6, zorder=4,
                        mfc=SURFACE if small else c, mec=c, mew=1.6)
    ax.set_yscale("log")
    ax.set_ylabel("Collapsed dwell (ns)", fontsize=9.5, color=INK)
    ax.set_xlabel("Sequence (designed k8 arrangements by number, see tools/tea_g6.py)",
                  fontsize=9, color=MUTED)
    _finish(fig, plt, "Collapsed-state dwells per sequence",
            rule + "\nOne faint point per dwell, jittered; solid marker = median, "
                   "bar = interquartile range.",
            [plt.Line2D([], [], marker="o", ls="", mfc=MUTED, mec=MUTED, ms=7,
                        label="median, bar = IQR"), open_marker],
            out_dir / "g6-random-dwells-strip.png")


# ---- main ------------------------------------------------------------------------

def _read(path: Path, text: tuple[str, ...]) -> list[dict]:
    return [{k: (v if k in text else float(v)) for k, v in r.items()}
            for r in csv.DictReader(open(path))]


def _write(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args(argv)
    run_csv = OUT_DIR / "g6-random.csv"
    ev_csv = OUT_DIR / "g6-random-dwells.csv"
    sum_csv = OUT_DIR / "g6-random-dwell-summary.csv"
    man_json = OUT_DIR / "g6-random-dwells.json"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        rows = _read(run_csv, ("run", "group", "label", "tea_gamma"))
        for r in rows:
            r["T"] = int(r["T"])
        events = _read(ev_csv, ("run", "rule"))
        summary = _read(sum_csv, ("run", "label", "rule"))
        for s in summary:
            s["T"] = int(s["T"])
        plot(rows, events, [s for s in summary if s["rule"] == "hyst"],
             json.loads(man_json.read_text()), OUT_DIR)
        return 0

    cells = [(g, lab, T) for g, (_, labs) in GROUPS.items() for lab in labs
             for T in TEMPERATURES]
    with Pool(args.jobs) as pool:
        results = pool.map(both, [run_name(lab, T) for _, lab, T in cells])
    missing = [run_name(lab, T) for (_, lab, T), r in zip(cells, results) if not r]
    if missing:
        print(f"not finished ({len(missing)}): {', '.join(missing)}")
    done = [(c, r) for c, r in zip(cells, results) if r]

    rng = np.random.default_rng(SEED)
    rows, events, summary = [], [], []
    for (g, lab, T), r in done:
        dw = r["dw"]
        for rule, mn in (("hyst", MIN_FRAMES), ("bare", 1)):
            d = np.array([L for L, _ in dw[rule]], dtype=float) * dw["dt_ns"]
            cens = np.array([c for _, c in dw[rule]], dtype=bool)
            events += [{"run": dw["run"], "rule": rule, "duration_ns": L, "censored": int(c)}
                       for L, c in zip(d, cens)]
            s = summarise(d, cens, mn, dw["dt_ns"], rng)
            summary.append({"run": dw["run"], "label": lab, "T": T, "rule": rule,
                            "collapse_fraction": dw["collapse_fraction"]} | s)
            if rule == "hyst":
                hyst = s
        rows.append({"group": g, "label": lab, "T": T} | r["rg"]
                    | {k: hyst[k] for k in ("n_dwells", "n_censored", "median_ns",
                                            "tau_ns", "tau_lo_ns", "tau_hi_ns")})

    _write(run_csv, rows)
    _write(ev_csv, events)
    _write(sum_csv, summary)
    dw0 = done[0][1]["dw"]
    dt = dw0["dt_ns"]
    manifest = {
        "runs_analysed": [r["run"] for r in rows], "runs_missing": missing,
        "groups": {g: labs for g, (_, labs) in GROUPS.items()},
        "n_residues": dw0["n_residues"], "rg_theta_nm": dw0["theta_nm"],
        "rg_theta_rule": "0.21 * N**0.5 nm",
        "analysis_window_ns": [DISCARD_NS, DISCARD_NS + dw0["frames"] * dt],
        "frame_ns": dt,
        "primary_rule": {"name": "hyst", "enter_below": f"{LOW} * Rg_theta",
                         "exit_above": f"{HIGH} * Rg_theta", "min_dwell_frames": MIN_FRAMES},
        "check_rule": {"name": "bare", "collapsed_if": "Rg < Rg_theta", "min_dwell_frames": 1},
        "min_dwell_ns": MIN_FRAMES * dt,
        "censoring": "dwells touching the first or last analysed frame are right-censored; "
                     "tau from censored geometric ML on frames above the minimum",
        "bootstrap": {"resamples": N_BOOT, "ci": "percentile 2.5-97.5", "seed": SEED},
        "small_n_flag": SMALL_N,
    }
    man_json.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {run_csv}, {ev_csv}, {sum_csv}, {man_json}")

    print(f"{'run':<24}{'Rg':>6}{'':>6}{'coll':>6}{'n':>5}{'tau':>7}{'tau 95% CI':>14}")
    for r in rows:
        print(f"{r['run'][3:]:<24}{r['rg_nm']:>6.2f} ±{r['rg_sem_nm']:.2f}"
              f"{r['collapse_fraction']:>6.0%}{r['n_dwells']:>5}{r['tau_ns']:>7.2f}"
              f"  {r['tau_lo_ns']:>5.2f}-{r['tau_hi_ns']:<5.2f}")
    plot(rows, events, [s for s in summary if s["rule"] == "hyst"], manifest, OUT_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
