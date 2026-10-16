#!/usr/bin/env bash
# The whole SubM4 analysis, end to end, from the repository root:   analysis/scripts/run_all.sh
# Step 01 extracts only runs that have no timeseries/<run>.npz yet; for many new runs submit
# analysis/scripts/extract.sbatch instead and rerun this script once it has finished.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
RUN=analysis/scripts/run.sh
$PY analysis/scripts/00_inventory.py
$PY analysis/scripts/01_extract_timeseries.py --workers "${WORKERS:-8}"
$RUN analysis/scripts/02_equilibration.py
$RUN analysis/scripts/03_metrics.py
$RUN analysis/scripts/04_stats.py
$RUN analysis/scripts/05_figures.py
$RUN analysis/scripts/06_reactive.py
$RUN analysis/scripts/07_models.py
$RUN analysis/scripts/05_figures.py --stats analysis/export/models/inter_log --out analysis/export/figures_logcoded --coding log
$RUN analysis/scripts/08_video_runs.py --workers "${WORKERS:-4}"
$RUN analysis/scripts/09_export.py
