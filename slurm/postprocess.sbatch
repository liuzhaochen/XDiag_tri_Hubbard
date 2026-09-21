#!/bin/bash
#SBATCH --job-name=triHub-post
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/postprocess-%j.out

set -euo pipefail

cd /home/zhl45qc/module/XDiag_tri_Hubbard
mkdir -p logs data/processed/N12 plots/N12
export MPLCONFIGDIR="${TMPDIR:-/tmp}/triHub-matplotlib-${SLURM_JOB_ID}"

julia --project=. --startup-file=no scripts/summarize_scan.jl \
    config/scan_coarse.toml data/processed/N12/coarse_summary.jld2
python3 scripts/plot_spectrum.py \
    data/processed/N12/coarse_summary.csv plots/N12
python3 scripts/analyze_spectrum.py \
    data/processed/N12/coarse_summary.csv data/processed/N12/coarse_analysis.txt
julia --project=. --startup-file=no scripts/summarize_scan.jl \
    config/scan_dense.toml data/processed/N12/dense_summary.jld2
python3 scripts/analyze_spectrum.py \
    data/processed/N12/dense_summary.csv data/processed/N12/dense_analysis.txt
python3 scripts/plot_crossing_zoom.py \
    data/processed/N12/dense_summary.csv plots/N12/crossing_zoom
