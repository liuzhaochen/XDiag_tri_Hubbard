# Implementation and calculation plan

## Development gates

- [x] Initialize the Julia project and persistent project ledger.
- [x] Validate the exact `3 x 4` triangular torus and all 36 bonds.
- [x] Generate triangular-lattice `D6` and recover the surviving torus point group.
- [x] Validate translations, `C2`, momentum stars, and little groups.
- [x] Construct XDiag `Electron(12,6,6,representation)` blocks.
- [x] Verify translation-sector dimensions reconstruct the raw dimension `853776`.
- [x] Verify `C2=+/-` dimensions reconstruct translation-only `(0,0)` and `(0,2)`.
- [x] Implement the Hubbard Hamiltonian and validate it against an independent
  fermionic bit-basis construction, including numerical symmetry commutators.
- [x] Reconstruct a small-cluster full spectrum from symmetry-resolved spectra.
- [x] Compare all five `k <-> -k` partner pairs at `U/t=6,8,12`.
- [x] Compare translation-only and `C2`-resolved spectra at `(0,0)` and `(0,2)`.
- [x] Benchmark multi-level Lanczos, reject its ghost spectrum, and validate guarded
  LOBPCG as the production solver.
- [ ] Run the 12-site coarse scan over `U/t = 6:0.5:14`.
- [ ] Analyze all sectors and the focused `Ky=0,pi` branches.
- [ ] Refine only around identified candidate crossings.
- [ ] Measure minimal observables only for states involved in a crossing.
- [ ] Discuss whether a `4 x 4` calculation is scientifically and computationally justified.

## Production policy

- The 12-site cluster is the only active cluster.
- Calculate seven inequivalent momentum stars, not all twelve momenta.
- Resolve `C2=+/-` only at `(0,0)` and `(0,2)`.
- Store compact JLD2 records; do not persist wavefunctions during the initial scan.
- One Slurm-array task owns one `(sector, U/t)` result.
- A result is accepted only after its solver diagnostics pass declared thresholds.
