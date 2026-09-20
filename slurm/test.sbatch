#!/bin/bash
#SBATCH --job-name=triHub-tests
#SBATCH --partition=test
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/test-%j.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs
julia --project=. --startup-file=no --check-bounds=yes test/runtests.jl
julia --project=. --startup-file=no scripts/report_sector_dimensions.jl
