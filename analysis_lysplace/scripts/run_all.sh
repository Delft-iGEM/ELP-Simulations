#!/usr/bin/env bash
# The RGD+ lysine-placement analysis end to end, from the repository root:  analysis_lysplace/scripts/run_all.sh
# Step 1 extracts only runs without timeseries/<run>.npz (~12 s each); for many new runs submit
# analysis_lysplace/scripts/extract.sbatch instead and rerun this script once it has finished.
set -euo pipefail
cd "$(dirname "$0")/../.."
RUN=analysis/scripts/run.sh
.venv/bin/python analysis_lysplace/scripts/00_inventory.py
.venv/bin/python analysis_lysplace/scripts/01_extract.py --workers "${WORKERS:-8}"
$RUN analysis_lysplace/scripts/02_metrics.py
$RUN analysis_lysplace/scripts/03_compare.py
$RUN analysis_lysplace/scripts/04_lysine_model.py
$RUN analysis_lysplace/scripts/05_figures.py
