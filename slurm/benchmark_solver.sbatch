#!/bin/bash
#SBATCH --job-name=triHub-solver
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/benchmark-solver-%j.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
julia --project=. --startup-file=no scripts/benchmark_multilevel_solver.jl
