#!/bin/bash
#SBATCH --job-name=sub4-extract
#SBATCH --partition=participation-cpu
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --output=analysis/logs/%x-%j.out
# Step 2a on a compute node: python analysis/scripts/01_extract_timeseries.py over all finished runs.
# Submit from the repository root:  sbatch analysis/scripts/extract.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1
srun "$SLURM_SUBMIT_DIR/.venv/bin/python" -u analysis/scripts/01_extract_timeseries.py --workers 16
