#!/bin/bash
#SBATCH --job-name=triHub-validate
#SBATCH --partition=test
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1-3
#SBATCH --output=logs/validate-sectors-%A-%a.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

u_values=(6.0 8.0 12.0)
u_value="${u_values[$((SLURM_ARRAY_TASK_ID - 1))]}"
julia --project=. --startup-file=no scripts/validate_sector_spectra.jl "${u_value}"
