# Spectral results

Updated: 2026-09-21

## Validated coarse scan

The scan contains six levels in all nine production blocks at
`U/t = 6:0.5:14` (153 independently validated records).

The absolute ground-state block on the coarse grid is:

- `(0,0,C2=+)` for `U/t=6.0,6.5`;
- `(1,0)` for `U/t=7.0,...,9.0`;
- `(0,0,C2=+)` for `U/t=9.5,...,14.0`.

Thus the two ground-branch crossings are bracketed by:

1. `6.5 < U/t < 7.0`;
2. `9.0 < U/t < 9.5`.

Both branches lie on the `Ky=0` line. The competing `(0,0)` state is even under
`C2`; the other branch has momentum `(1,0)` and represents the star
`{(1,0),(2,0)}`. The `Ky=pi` branches do not become the ground state on this grid.

The smallest coarse-grid many-body gap is `0.00225145631249 t` at `U/t=6.5`,
consistent with the first crossing lying close to that grid point.

Several additional crossings occur between excited block-lowest branches. They are
not refinement targets until the two ground crossings are resolved.

## Completed dense refinement

A spacing-0.05 scan is restricted to `(0,0,C2=+)` and `(1,0)` at
`U/t=6.55:0.05:6.95` and `9.05:0.05:9.45` (36 records total).
All records passed the numerical gates. The two sampled ground-crossing brackets
are:

- `6.60 < U/t < 6.65`, with linear energy-difference zero near `6.6281`;
- `9.30 < U/t < 9.35`, with linear energy-difference zero near `9.3124`.

The interpolated values are visual guides, not additional ED resolution and not
thermodynamic transition estimates.

An all-12-momentum translation-only audit independently confirmed both order
reversals and the exact `(1,0)<->(2,0)` partner degeneracy. At `U/t=8`, the pair
lies `0.033449734171 t` below the complete `(0,0)` block. A separate `Sz=1` scan
lies `0.077558948956 t` above it, so the finite-momentum ground is a singlet.

## Other crossing candidates

- The lowest `(0,0,C2=-)` and `(1,0)` branches cross in the coarse brackets
  `6.0--6.5` and `13.0--13.5`; the first is very low lying but not the ground.
- On `Ky=pi`, the lowest `(0,2,C2=-)` and `(1,2)` branches cross in
  `9.5--10.0`, only about `0.01--0.02 t` above the ground at the endpoints.
- Numerous higher retained levels exchange order across distinct irreps. They are
  crossing candidates, not evidence of thermodynamic phases.

The reader-facing interpretation is in `docs/physical_report.md`.
