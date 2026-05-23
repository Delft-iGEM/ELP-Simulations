#!/bin/bash
#SBATCH --job-name={sim_name}
#SBATCH --partition=gpu-a100
#SBATCH --dependency=singleton
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --output=logs/%x-%j.out

set -euo pipefail
module purge
module load 2024r1
module load cuda

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

srun uv run python -u -m "simulations.{sim_name}.runtime.run"
