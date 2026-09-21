# Implementation and calculation plan

## Development gates

- [x] Initialize the Julia project and persistent project ledger.
- [x] Validate the exact `3 x 4` triangular torus and all 36 bonds.
- [x] Generate triangular-lattice `D6` and recover the surviving torus point group.
- [x] Validate translations, `C2`, momentum stars, and little groups.
- [x] Construct XDiag `Electron(12,6,6,representation)` blocks.
- [x] Verify translation-sector dimensions reconstruct the raw dimension `853776`.
- [x] Verify the `A/B` little-group irreps reconstruct translation-only `Gamma`
  and `M2` blocks.
- [x] Implement the Hubbard Hamiltonian and validate it against an independent
  fermionic bit-basis construction, including numerical symmetry commutators.
- [x] Reconstruct a small-cluster full spectrum from symmetry-resolved spectra.
- [x] Compare all five `k <-> -k` partner pairs at `U/t=6,8,12`.
- [x] Compare translation-only and `C2`-resolved spectra at `(0,0)` and `(0,2)`.
- [x] Benchmark multi-level Lanczos, reject its ghost spectrum, and validate guarded
  LOBPCG as the production solver.
- [x] Run and validate the 12-site coarse scan over `U/t = 6:0.5:14`.
- [x] Analyze all sectors and the focused `Ky=0,pi` branches.
- [x] Refine the two identified ground crossings on the `Ky=0` line.
- [x] Rerun all 153 coarse-grid blocks and measure double occupancy, spin/charge
  correlations, and uniform chirality for all six retained levels.
- [x] Measure the ground-state scalar-chirality structure factor at all 12 allowed
  momenta, retaining the full up/down-triangle sublattice Gram matrix.
- [x] Generate ground-state, low-state, representative heatmap, and cluster
  schematic plots; incorporate their interpretation in the physical report.
- [x] Estimate a possible `4 x 4` calculation and defer it until explicitly
  requested. A future implementation must use the full surviving `D6` symmetry,
  not reuse the `3 x 4` `C2` sector inventory.

## Production policy

- The 12-site cluster is the only active cluster.
- Calculate seven inequivalent momentum stars, not all twelve momenta.
- Resolve the `C2` little-group irreps `A/B` only at `Gamma=(0,0)` and
  `M2=(0,2)`.
- Store compact JLD2 records; do not persist wavefunctions during the initial scan.
- One Slurm-array task owns one `(sector, U/t)` result.
- A result is accepted only after its solver diagnostics pass declared thresholds.
- Continuous heatmap colors are finite Fourier visualizations; allowed cluster
  momenta are always marked explicitly as dots.
- A `4 x 4` study is not active production work. If resumed, begin with a
  full-`D6` sector implementation and a one-point resource benchmark before any
  17-point scan.
