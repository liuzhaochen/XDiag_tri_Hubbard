# Validation ledger

Updated: 2026-09-21

## Required geometry and symmetry checks

- [x] 12 sites with the specified site ordering.
- [x] 36 unique undirected nearest-neighbor bonds.
- [x] Every site has coordination number 6.
- [x] Triangular `D6` contains exactly 12 integer-coordinate operations.
- [x] The `3 x 4` torus retains exactly identity and inversion `C2`.
- [x] Every retained operation preserves the bond set.
- [x] The translation and point-permutation groups close with `Tx^3=Ty^4=C2^2=1`.
- [x] The 12 momenta form exactly seven `C2` stars.

## Required XDiag checks

- [x] Julia 1.13.0 and XDiag 0.5.1 recorded from a compute node.
- [x] All 12 translation-block dimensions sum to `binomial(12,6)^2 = 853776`.
- [x] Star-weighted dimensions of the seven representatives give the same total.
- [x] `C2=+/-` dimensions sum to the translation-only dimension at `(0,0)` and `(0,2)`.
- [x] A complete six-site Hubbard spectrum agrees with the union of all translation
  blocks at `U/t=3.7` to absolute tolerance `1e-11`.
- [x] All five partner pairs have matching six-level spectra at `U/t=6,8,12`.
- [x] Multilevel Lanczos was explicitly rejected; guarded LOBPCG passes residual and
  orthogonality gates on the largest block at `U/t=8`.

The compute-node suite passes 46/46 tests. The direct `HubbardU` convention
was independently checked on the complete two-site `(Nup,Ndn)=(1,1)` basis.
An independent fermionic bit-basis implementation verifies numerical commutators
with `Tx`, `Ty`, and `C2` on a six-site torus to `1e-13`, and its complete Hubbard
spectrum agrees with XDiag to `1e-11`.

## Measured multilevel solver benchmark

- Sector `(1,0)`, dimension 71,178, `U/t=8`, four CPUs.
- CSR32: 2,866,586 nonzeros, 57.6 MB Julia summary size, built in 1.72 s.
- Six-state Lanczos: 197 iterations and 8.15 s, but maximum explicit residual
  `1.11e-1` and orthogonality error `9.996e-1`; rejected despite criterion
  `converged`.
- Six-state LOBPCG plus two guards: 284 iterations and 33.51 s; maximum reported
  residual `9.84e-11`, explicit residual agreement at roundoff, orthogonality error
  `1.78e-15`.
- The trustworthy ground values agree to `3.6e-14`.

## Completed production-sector validation

- Slurm array job 3324639 completed all three tasks at `U/t=6,8,12`.
- All five `k <-> -k` pairs agree for six low-energy levels; the largest absolute
  discrepancy is `1.95e-14 t`.
- At `(0,0)` and `(0,2)`, the first six translation-only levels agree with the
  sorted union of the `C2=+/-` levels at all three couplings, with the same
  `1.95e-14 t` worst-case scale.
- Task elapsed times were 8:43--14:10 and peak RSS was 588,664 KiB.

## Completed coarse scan

- Slurm array job 3325184 produced 152 accepted records directly. Task 153
  `(1,3), U/t=14` exhausted three 500-iteration seeds and wrote no partial result.
- A targeted 1000-iteration rerun converged on its first seed in 74.24 s with
  maximum residual `9.66e-11`, completing all 153 expected records.
- The summary validator found exactly 153 records with no unexpected files.
- Maximum reported residual: `8.59e-9`; maximum independently recomputed residual:
  `8.59e-9`; maximum eigenvector orthogonality error: `3.00e-15`.
- Array task walltimes were 0:27--4:04 and maximum peak RSS was 589,116 KiB.
- Portable outputs are `data/processed/N12/coarse_summary.jld2` and
  `coarse_summary.csv`; five overview plots were generated as PDF and PNG.

## Completed dense crossing scan

- Slurm array job 3325410 completed all 36 targeted records in the
  `(0,0,C2=+)` and `(1,0)` blocks.
- The validator found exactly 36 records with no unexpected files.
- Maximum reported residual: `1.49e-10`; maximum independently recomputed residual:
  `1.49e-10`; maximum eigenvector orthogonality error: `2.89e-15`.
- Array task walltimes were 0:26--1:38 and maximum peak RSS was 557,180 KiB.
- The two `0.05` crossing brackets are `6.60--6.65` and `9.30--9.35`.
