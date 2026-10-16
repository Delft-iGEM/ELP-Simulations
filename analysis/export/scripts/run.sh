#!/usr/bin/env bash
# Run an analysis script with the project venv plus the analysis-only packages, without
# touching pyproject.toml / uv.lock:   analysis/scripts/run.sh analysis/scripts/04_stats.py
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"
exec uv run --no-sync --with scikit-learn==1.5.2 --with plotly==5.24.1 --with kaleido==0.2.1 \
    --with shap==0.46.0 python "$@"
