"""Step 9 - assemble analysis/export/ for the wiki write-up: copies, MANIFEST.md, README.md, zip.

Run after 07_models.py (models/, needs 05_figures.py --coding log for figures_logcoded/) and
08_video_runs.py (video_runs*). Copies the finished analysis outputs next to them, writes
MANIFEST.md (every file, what it holds, its columns, the date and the commit / hashes of the
scripts) and the handoff README.md, then analysis_export.zip with every file of at most 25 MB.
No trajectories or .npz files are exported.

    analysis/scripts/run.sh analysis/scripts/09_export.py
"""
from __future__ import annotations

import datetime as dt
import fnmatch
import gzip
import hashlib
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from common import ANALYSIS, ROOT  # noqa: E402

OUT = ANALYSIS / "export"
ZIP = OUT / "analysis_export.zip"
MAX_BYTES = 25 * 1024 ** 2
COPIES = {                                    # source (relative to analysis/) -> export path
    "stats": "stats", "figures": "figures", "scripts": "scripts",
    "runs_metrics.csv": "runs_metrics.csv", "run_inventory.csv": "run_inventory.csv",
    "per_lysine.csv": "per_lysine.csv", "timeseries_primary.csv.gz": "timeseries_primary.csv.gz",
    "reactive_anchors.csv": "reactive_anchors.csv", "settings.yaml": "settings.yaml",
    "qc/equilibration_per_run.csv": "qc/equilibration_per_run.csv", "README.md": "analysis_README.md",
}
DESCRIBE = [                                  # first matching pattern wins
    ("README.md", "Handoff: what is here, what is new, what could not be done, how to regenerate."),
    ("MANIFEST.md", "This file."),
    ("analysis_README.md", "The full analysis write-up (analysis/README.md): methods, QC, results, open questions."),
    ("settings.yaml", "Frozen definitions (response, pair definition, cutoff, burn-in, QC thresholds, coding ranges, seed)."),
    ("runs_metrics.csv", "One row per run: design parameters, all per-run metrics with block SEs, QC flags. "
                         "`data_source` = trajectory (104) or teammate_summary (20, team-definition-A columns only)."),
    ("run_inventory.csv", "Every run-sheet row matched to its output, with status (finished / summary_only / excluded / "
                          "missing) and the reason."),
    ("per_lysine.csv", "Per lysine per run: chain, residue (0-based in chain), K index, pinned, mean free inter-chain partners at 12 Å."),
    ("timeseries_primary.csv.gz", "Per-frame inter / intra free contacts at 12 Å for every trajectory run, all frames including burn-in."),
    ("reactive_anchors.csv", "Titus's reactive anchors A10-A19: team-A contacts and crosslinks formed / loops / bridges. Never in the contact fits."),
    ("qc/equilibration_per_run.csv", "Chodera burn-in per run and observable (frame index)."),
    ("stats/summary.txt", "Frozen first-pass statistics in plain text (n = 104)."),
    ("stats/coefficients.csv", "Frozen model: OLS coefficients on log(inter_free_12A), SE, t, p, 95 % CI, VIF."),
    ("stats/parameter_tests.csv", "Frozen model: per-parameter F-tests and LMG shares of R²."),
    ("stats/lack_of_fit.csv", "Frozen model: lack-of-fit F-test against the anchor pure error."),
    ("stats/model_fit.csv", "Frozen model: R², grouped CV (OLS and random forest), residual tests, Cook's distance."),
    ("stats/run_diagnostics.csv", "Frozen model, per run: fitted, residual, leverage, Cook's D, held-out predictions."),
    ("stats/importance_rf.csv", "Random-forest permutation importance and mean |SHAP| (frozen analysis)."),
    ("stats/shap_values.npy", "SHAP values of the random forest (NumPy array, runs x parameters)."),
    ("stats/glm_gamma.csv", "Gamma GLM (log link) coefficients next to OLS on log y."),
    ("stats/sensitivity.csv", "The frozen model on the sensitivity responses (8 / 16 Å, team definition A with and without Titus's runs, 200 ns burn-in)."),
    ("stats/secondary.csv", "The frozen model on secondary metrics (R² and LMG shares)."),
    ("stats/teammate_check.csv", "Titus's 20 runs predicted by the team-A model fitted on our runs alone (held out)."),
    ("stats/reactive_anchors.csv", "Reactive anchors vs the non-reactive team-A model at the same design points."),
    ("stats/reactive_glm.csv", "Reactive anchors: quasi-binomial GLM of the crosslink bridge fraction (main effects)."),
    ("stats/reactive_design_points.csv", "Reactive anchors per design point, with non-reactive model predictions."),
    ("stats/reactive_summary.txt", "Reactive-anchor analysis in plain text."),
    ("models/comparison.csv", "All six new fits side by side: R², adjusted R², grouped CV, lack of fit, residual tests, Cook's distance."),
    ("models/lmg_shares.csv", "Per-parameter F-tests and LMG shares for all six fits."),
    ("models/*/coefficients.csv", "Coefficients, SE, t (OLS) or z (GLM), p, 95 % CI, VIF, fold change / odds ratio per coded unit."),
    ("models/*/parameter_tests.csv", "Per-parameter F-test (all terms containing the parameter dropped) and LMG shares."),
    ("models/*/lack_of_fit.csv", "Lack-of-fit F-test against the triplicate anchors' pure error (log scale, or deviance for the GLM)."),
    ("models/*/model_fit.csv", "R², adjusted R² (deviance-based for the GLM), 5-fold grouped CV, residual tests, VIF, Cook's distance."),
    ("models/*/run_diagnostics.csv", "Per run: observed, fitted (link scale), residual, studentised residual, leverage, Cook's D, held-out prediction."),
    ("models/*/effects_oat.csv", "One parameter at a time, others at the base case: low / base / high predictions with 95 % CI and the high/low ratio."),
    ("models/*/importance_rf.csv", "Random-forest importances copied from stats/ (trees are invariant to a monotone recoding)."),
    ("figures/*", "Frozen-model figures (Plotly HTML with hover, and PNG); captions in analysis_README.md §8."),
    ("figures_logcoded/*", "The same figures for the log-coded primary model."),
    ("video_runs.md", "The runs chosen for the videos: reasons, absolute paths, production start frame, frame count."),
    ("video_runs/*_contacts.csv", "Every frozen-definition contact per frame (12 Å, free lysines, all frames); see video_runs.md."),
    ("video_runs/*_pins.csv", "Surface-pinned lysines: chain, resid, K index, atom index."),
    ("video_runs/*_top.pdb", "Copy of the run's topology (Å)."),
    ("video_runs/*_surface_state.json", "Copy of the run's surface state (pinned lysines, tether points)."),
    ("scripts/*", "Copy of the analysis scripts that produced everything here."),
]


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def copy_outputs():
    for src, dst in COPIES.items():
        s, d = ANALYSIS / src, OUT / dst
        if s.is_dir():
            shutil.rmtree(d, ignore_errors=True)
            shutil.copytree(s, d, ignore=shutil.ignore_patterns("__pycache__"))
        else:
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(s, d)


def columns(p: Path) -> str:
    if p.name.endswith(".csv.gz"):
        with gzip.open(p, "rt") as f:
            head = f.readline()
    elif p.suffix == ".csv":
        with open(p, encoding="utf-8-sig") as f:
            head = f.readline()
    else:
        return ""
    return ", ".join(f"`{c}`" for c in head.strip().split(","))


def describe(rel: str) -> str:
    for pat, text in DESCRIBE:
        if fnmatch.fnmatch(rel, pat):
            return text
    return ""


def files():
    return sorted(p for p in OUT.rglob("*") if p.is_file() and p != ZIP)


def write_readme(now, head, dirty):
    comp = pd.read_csv(OUT / "models" / "comparison.csv").set_index(["response", "coding"])
    lin, log = comp.loc[("inter", "linear")], comp.loc[("inter", "log")]
    oat = pd.read_csv(OUT / "models" / "inter_log" / "effects_oat.csv").set_index("parameter")
    oatl = pd.read_csv(OUT / "models" / "inter_linear" / "effects_oat.csv").set_index("parameter")
    inv = pd.read_csv(ANALYSIS / "run_inventory.csv")
    diag = pd.read_csv(OUT / "models" / "inter_log" / "run_diagnostics.csv")
    top = diag.loc[diag.cooks_d.idxmax()]

    def row(r, k):
        return (f"| {k[0]} | {k[1]} | {r.n_runs} | {r.R2:.3f} | {r.adj_R2:.3f} | {r.cv5_R2:.3f} | "
                f"{r.cv5_RMSE_link:.3f} | {r.lof_F:.2f} ({r.lof_df_num:.0f}, {r.lof_df_den:.0f}) | {r.lof_p:.2g} |")

    sig = "still significant" if log.lof_p < 0.05 else "no longer significant"
    lines = [
        "# SubM4 export for the write-up", "",
        f"Generated {now} by `analysis/scripts/09_export.py`, repository HEAD `{head}`. "
        + ("**The analysis scripts are not committed** (`analysis/` is untracked), so `MANIFEST.md` "
           "lists a SHA-256 for each script and `scripts/` holds a copy." if dirty else ""), "",
        "## What is here", "",
        "* The finished SubM4 analysis: `stats/`, `figures/`, `runs_metrics.csv`, `run_inventory.csv`, "
        "`per_lysine.csv`, `timeseries_primary.csv.gz`, `reactive_anchors.csv`, `settings.yaml`, "
        "`qc/equilibration_per_run.csv`, and the full write-up as `analysis_README.md`.",
        "* `models/`: six response-surface fits, each with coefficients and 95 % CIs, fit statistics, grouped CV, "
        "lack of fit, residual diagnostics, Cook's distance, the one-parameter-at-a-time table and LMG shares "
        "(`comparison.csv` and `lmg_shares.csv` put them side by side). Responses: `inter` = log inter-chain "
        "free contacts (frozen primary), `intra` = log intra-chain free contacts, `interfrac` = inter / (inter + "
        "intra) as a quasi-binomial GLM with logit link on the production pair-frame counts. Each has a `linear` "
        "(frozen) and a `log` coding. `models/inter_linear` reproduces `stats/` exactly.",
        "* `figures_logcoded/`: the frozen figure set redrawn for the log-coded primary model.",
        "* `video_runs.md` and `video_runs/`: the runs picked for the videos, with paths, the production start "
        "frame, and per-frame contact and pin files (checked frame by frame against the metrics).",
        "* `MANIFEST.md`: every file with its contents and columns. `analysis_export.zip`: all of it.", "",
        "## New: the log-coded model", "",
        "The full quadratic was refitted with the natural log of all four parameters, each coded to −1…+1 over "
        f"its log range. On the primary response it raises R² from {lin.R2:.3f} to {log.R2:.3f} and grouped-CV R² "
        f"from {lin.cv5_R2:.3f} to {log.cv5_R2:.3f}, cuts the lack of fit from F = {lin.lof_F:.2f} to "
        f"{log.lof_F:.2f} ({sig}, p = {log.lof_p:.2g}), and removes the spurious density maximum: the "
        f"density effect is ×{oat.loc['concentration', 'ratio_high_over_low']:.1f} "
        f"[{oat.loc['concentration', 'ratio_ci_low']:.1f}, {oat.loc['concentration', 'ratio_ci_high']:.1f}] "
        f"from 0.01 to 0.06 chains/nm², against ×{oatl.loc['concentration', 'ratio_high_over_low']:.1f} in the frozen "
        "model, and it rises monotonically.", "",
        "**The log coding was chosen after the lack of fit of the frozen model had been seen.** It is not part of "
        "the frozen analysis plan. The frozen linear-coded model stays the baseline, and any claim based on the log "
        "coding should say so.", "",
        f"The most influential run under the log coding is `{top.slug}` (Cook's D = {top.cooks_d:.2f}, studentised "
        f"residual {top.resid_studentized:.1f}), the low-density, narrow-spacing corner already flagged in the "
        "frozen fit; it passes QC and is kept.", "",
        "| response | coding | n | R² | adj. R² | CV R² | CV RMSE (link) | lack of fit F (df) | p |",
        "|---|---|---|---|---|---|---|---|---|",
        *[row(comp.loc[k], k) for k in comp.index], "",
        "R² for `interfrac` is the deviance pseudo-R²; its lack of fit is deviance-based and its residual tests use "
        "deviance (Shapiro–Wilk) and Pearson (Breusch–Pagan) residuals. The dispersion is estimated, so the "
        "autocorrelated frames are not treated as independent trials.", "",
        "## What could not be done, and choices made", "",
        "* **Wet-lab constructs (video item 6): none selected.** No Monoblock, Triblock-32, Triblock-64 or "
        "Tetrablock run is in the SubM4 dataset. Folders with those names exist from earlier submodules, but "
        "none has a readable trajectory here (see `video_runs.md` §6).",
        "* **Git commit of the scripts:** there is none, because `analysis/` is untracked. Hashes are in `MANIFEST.md`.",
        f"* **Titus's {int((inv.status == 'summary_only').sum())} `sub4-extra` runs** have summaries only, not "
        "trajectories, so they are in none of these fits (they lack the free-only and intra-chain counts). They "
        "enter only the team-definition-A sensitivity fit in `stats/`.",
        "* **\"Bridge fraction\"** is taken as the non-reactive inter / (inter + intra) free-contact fraction "
        "(`interfrac`). The crosslink bridge fraction of the reactive anchors A10–A19 has 10 design points, too few "
        "for a 15-term quadratic; its main-effects analysis is in `stats/reactive_*`.",
        "* **Random forest** importances in `figures_logcoded/importance_comparison` are the frozen ones, not a "
        "refit: trees split on order, so a monotone recoding of the inputs gives the same partitions.",
        "* **Video item 4:** no run has both the longest chain and the widest spacing; the choice rule and the "
        "two single-criterion runs are in `video_runs.md`.",
    ]
    big = [p for p in files() if p.stat().st_size > MAX_BYTES]
    lines.append("* **Zip:** every file is under 25 MB, nothing left out." if not big else
                 "* **Left out of the zip (> 25 MB):** " + ", ".join(f"`{p.relative_to(OUT)}`" for p in big))
    pending = sorted(p.name for p in (ROOT / "simulations").glob("sub4-*")
                     if p.is_dir() and not (inv.slug == p.name).any())
    if pending:
        lines.append(f"* **Not yet included:** {len(pending)} prepared SubM4 runs that are in no sheet listed in "
                     "`scripts/common.py` (e.g. the `sub4-main-*` runs of `Chr_sub4_run5.csv` and the `-nr` A10–A19 "
                     "reruns). Add their sheets to `LOCAL_SHEETS` when they finish and rerun the commands below.")
    lines += ["", "## Regenerate", "",
              "From the repository root on DelftBlue (`run.sh` layers the pinned analysis packages on the project venv):",
              "", "```bash", "analysis/scripts/run_all.sh            # steps 00-09 in order", "```", "",
              "or step by step:", "", "```bash",
              ".venv/bin/python analysis/scripts/00_inventory.py",
              ".venv/bin/python analysis/scripts/01_extract_timeseries.py --workers 8   # only runs without a .npz; or sbatch analysis/scripts/extract.sbatch",
              "analysis/scripts/run.sh analysis/scripts/02_equilibration.py",
              "analysis/scripts/run.sh analysis/scripts/03_metrics.py",
              "analysis/scripts/run.sh analysis/scripts/04_stats.py",
              "analysis/scripts/run.sh analysis/scripts/05_figures.py",
              "analysis/scripts/run.sh analysis/scripts/06_reactive.py",
              "analysis/scripts/run.sh analysis/scripts/07_models.py",
              "analysis/scripts/run.sh analysis/scripts/05_figures.py --stats analysis/export/models/inter_log --out analysis/export/figures_logcoded --coding log",
              "analysis/scripts/run.sh analysis/scripts/08_video_runs.py",
              "analysis/scripts/run.sh analysis/scripts/09_export.py", "```", "",
              "New runs: add the sheet to `LOCAL_SHEETS` in `scripts/common.py`, then run everything. "
              "Nothing in `settings.yaml` changes: the burn-in and every other definition stay frozen. "
              "The video selection is re-derived by rule, so it can change when new runs arrive (e.g. a new lowest density)."]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def write_manifest(now, head, dirty, scripts):
    lines = ["# Manifest — SubM4 export", "",
             f"Generated {now} by `analysis/scripts/09_export.py`. Repository HEAD: `{head}` "
             f"({git('log', '-1', '--format=%ci', 'HEAD')}).", "",
             ("The analysis scripts are **not committed** (`analysis/` is untracked in git), so no commit identifies "
              "them; their SHA-256 hashes at export time are listed below, and `scripts/` holds the copies." if dirty
              else "The analysis scripts are committed at the HEAD above."), "",
             "| script | SHA-256 |", "|---|---|",
             *[f"| `{n}` | `{h}` |" for n, h in scripts], "",
             "## Files", "",
             "Paths are relative to `export/`. Sizes in kB. Columns are given for CSV files.", ""]
    for p in files():
        rel = str(p.relative_to(OUT))
        lines.append(f"### `{rel}` ({p.stat().st_size / 1024:.0f} kB)")
        lines.append("")
        d = describe(rel)
        if d:
            lines.append(d)
        c = columns(p)
        if c:
            lines += ["", f"Columns: {c}"]
        lines.append("")
    lines += ["### `analysis_export.zip`", "", "Everything above in one archive (files over 25 MB left out; see README.md).", ""]
    (OUT / "MANIFEST.md").write_text("\n".join(lines))


def write_zip():
    skipped = []
    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in files():
            if p.stat().st_size > MAX_BYTES:
                skipped.append(p)
                continue
            z.write(p, Path("analysis_export") / p.relative_to(OUT))
    return skipped


def main():
    for need in ["models/comparison.csv", "figures_logcoded/coefficients.html", "video_runs.md"]:
        if not (OUT / need).exists():
            sys.exit(f"export/{need} missing - run 07_models.py, 05_figures.py --coding log and 08_video_runs.py first")
    copy_outputs()
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    head = git("rev-parse", "--short", "HEAD")
    dirty = bool(git("status", "--porcelain", "--", "analysis/scripts"))
    scripts = [(p.name, hashlib.sha256(p.read_bytes()).hexdigest())
               for p in sorted((ANALYSIS / "scripts").iterdir()) if p.is_file()]
    (OUT / "MANIFEST.md").touch()           # so the file list includes it on a first run too
    write_readme(now, head, dirty)
    write_manifest(now, head, dirty, scripts)
    skipped = write_zip()
    n = len(files())
    print(f"{n} files in export/; zip {ZIP.stat().st_size / 1024 ** 2:.1f} MB; "
          f"left out: {[str(p.relative_to(OUT)) for p in skipped] or 'nothing'}")


if __name__ == "__main__":
    main()
