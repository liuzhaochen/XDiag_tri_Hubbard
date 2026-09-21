#!/bin/bash
#SBATCH --job-name=triHub-dense
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1-36%36
#SBATCH --output=logs/dense-%A-%a.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

sector_position=$(( (SLURM_ARRAY_TASK_ID - 1) % 2 + 1 ))
parameter_index=$(( (SLURM_ARRAY_TASK_ID - 1) / 2 + 1 ))
julia --project=. --startup-file=no scripts/run_scan.jl \
    config/scan_dense.toml "${sector_position}" "${parameter_index}"
