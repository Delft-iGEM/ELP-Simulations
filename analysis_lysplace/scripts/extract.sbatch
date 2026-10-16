#!/bin/bash
#SBATCH --job-name=lysplace-extract
#SBATCH --partition=participation-cpu
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --output=analysis_lysplace/logs/%x-%j.out
# Step 1 on a compute node. Submit from the repository root:  sbatch analysis_lysplace/scripts/extract.sbatch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1
srun "$SLURM_SUBMIT_DIR/.venv/bin/python" -u analysis_lysplace/scripts/01_extract.py --workers 16
