# XDiagTriHubbard

Symmetry-resolved exact diagonalization of the half-filled `3 x 4`
triangular-lattice Hubbard model with XDiag.

The initial scope is the 12-site torus with `Nup = Ndn = 6`. A `4 x 4` cluster is
explicitly deferred until the complete 12-site spectrum and measured resource costs
have been reviewed.

Julia must never run on the login node. Run tests on a Slurm compute node:

```bash
sbatch slurm/test.sbatch
```

or interactively:

```bash
srun --ntasks=1 --immediate=60 -p test --time=00:30:00 --cpus-per-task=1 \
  julia --project=. --startup-file=no --check-bounds=yes test/runtests.jl
```

Persistent project status, decisions, and validation evidence live in
[`.agents/`](.agents/README.md). The scientific specification is
[`.agents/prompt.md`](.agents/prompt.md).

The finite-cluster setup, symmetry blocks, numerical evidence, observed crossings,
and plot index are collected in the
[`physical report`](docs/physical_report.md).

The coarse production grid is declared in
[`config/scan_coarse.toml`](config/scan_coarse.toml). Each Slurm-array task writes
one atomic raw record:

```bash
sbatch slurm/coarse_scan.sbatch
```

After all 153 records pass their local solver gates, validate and aggregate them:

```bash
srun --ntasks=1 --immediate=60 -p test --time=00:15:00 --cpus-per-task=1 \
  julia --project=. --startup-file=no scripts/summarize_scan.jl \
  config/scan_coarse.toml data/processed/N12/coarse_summary.jld2
```

The portable raw and processed schemas are documented in
[`docs/data_format.md`](docs/data_format.md).

The validated coarse-grid observable workflow uses two arrays: one for all six
retained states in every block and one for the ground-state two-sublattice
chirality structure factor:

```bash
sbatch slurm/observable_scan.sbatch
sbatch slurm/ground_chirality.sbatch
```

Aggregate only after all 170 raw records exist:

```bash
srun --ntasks=1 --immediate=60 -p test --time=00:15:00 --cpus-per-task=1 \
  julia --project=. --startup-file=no scripts/summarize_observables.jl \
  config/scan_observables.toml data/processed/N12/observables_summary.jld2
```

Plot commands and the distinction between exact allowed-momentum data and
continuous Fourier visualization are documented in [`plots/README.md`](plots/README.md).
