#!/bin/bash
#SBATCH --job-name=triHub-obs
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1-153%36
#SBATCH --output=logs/observables-%A-%a.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

sector_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 9 + 1 ))
parameter_index=$(( (SLURM_ARRAY_TASK_ID - 1) / 9 + 1 ))
julia -t "${SLURM_CPUS_PER_TASK}" --project=. --startup-file=no \
    scripts/run_observables.jl config/scan_observables.toml \
    "${sector_index}" "${parameter_index}"
