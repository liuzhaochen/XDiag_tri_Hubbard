#!/bin/bash
#SBATCH --job-name=triHub-chi
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1-17%17
#SBATCH --output=logs/chirality-%A-%a.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
julia -t "${SLURM_CPUS_PER_TASK}" --project=. --startup-file=no \
    scripts/run_ground_chirality.jl config/scan_observables.toml \
    "${SLURM_ARRAY_TASK_ID}"
