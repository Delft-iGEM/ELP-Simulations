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

module load 2026
module load cuda/12.9

export PATH="$HOME/.local/bin:$PATH"

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

# Call the venv's interpreter directly, NOT `uv run`. `uv run` re-syncs the
# environment on every invocation, and when a batch of jobs starts at once one
# can read site-packages while another is mid-rewrite. That cost three runs out
# of 125 launches (2.4%), each dying within seconds on a missing matplotlib data
# file, e.g.
#     FileNotFoundError: .../mpl-data/stylelib/seaborn-v0_8-paper.mplstyle
# The file is always there afterwards, which is what gives the race away. The
# environment is already synced by `sim new` / `sim prepare` on the login node,
# so there is nothing for the compute node to sync.
srun "$SLURM_SUBMIT_DIR/.venv/bin/python" -u -m "simulations.{sim_name}.runtime.run"
