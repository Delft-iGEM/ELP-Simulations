#!/bin/bash
#SBATCH --job-name={sim_name}
#SBATCH --partition={partition}
#SBATCH --dependency=singleton
#SBATCH --time={runtime}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpu_per_task}
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=8000
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

module load 2024r1
module load cuda

export PATH="$HOME/.local/bin:$PATH"

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

srun uv run python -u -m "simulations.{sim_name}.runtime.run"
