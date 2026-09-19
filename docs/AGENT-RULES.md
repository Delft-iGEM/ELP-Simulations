# Rules for anyone (human or AI) working on the robust-audit branch

* Work ONLY inside /home/tnartey/ELP-Simulations-robust (git worktree, branch robust-audit).
  NEVER edit /home/tnartey/ELP-Simulations (main): 20 SLURM jobs import tools/ from it.
* Python: /home/tnartey/ELP-Simulations/.venv/bin/python, run with cwd = the worktree
  (cwd is first on sys.path so the worktree's tools/ wins). Never `uv run` here (it would
  build a second venv). Never `sbatch`. CPU only, <= 4 threads, keep test systems tiny.
* Test data: export ELP_DATA_DIR=/scratch/tnartey/elp-audit-data (runtime/ symlinks go there).
* Test simulations are named zz-audit-*; keep them small (<= 4 chains, <= 60 residues,
  <= 20000 steps) — this is a login node.
* Every bug found goes in docs/BUGS.md (id, file:line, symptom, root cause, fix/status).
* Every modelling gap or opportunity goes in docs/GAPS.md.
