# Finite-cluster level-crossing report

Updated: 2026-09-21

## Purpose and interpretation

This calculation searches for low-energy level crossings in one small,
symmetry-resolved Hubbard cluster. It is **not** an attempt to reproduce an iDMRG
phase diagram with a 12-site system. In particular, the crossing locations below
are finite-cluster spectral diagnostics, not thermodynamic transition estimates.

## Model and cluster

The Hamiltonian is

\[
H=-t\sum_{\langle ij\rangle,\sigma}
(c^\dagger_{i\sigma}c_{j\sigma}+\mathrm{h.c.})
+U\sum_i n_{i\uparrow}n_{i\downarrow},\qquad t=1.
\]

We use a half-filled `3 x 4` triangular torus with 12 sites and
`Nup=Ndn=6`. Its periods are `3 a1` and `4 a2`, and the three positive bond
directions are `a1`, `a2`, and `a2-a1`. The generated graph has 36 undirected
bonds and coordination six. Before symmetry reduction, the fixed-particle-number
Hilbert space has

\[
\binom{12}{6}^2=853,776
\]

states.

## Symmetry and block labels

Translations obey `Tx^3=Ty^4=1`, so momenta are labeled by
`(nx,ny)`, with `nx=0,1,2`, `ny=0,1,2,3`, and
`Ky=2 pi ny/4`. Mechanical reduction of triangular-lattice `D6` on this torus
leaves only

\[
C_2:(x,y)\mapsto(-x,-y),\qquad
(n_x,n_y)\mapsto(-n_x,-n_y).
\]

Thus the 12 translation momenta form seven `k <-> -k` stars. We diagonalize one
representative of each star. Only the self-conjugate momenta `(0,0)` and `(0,2)`
have an additional one-dimensional little-group label `C2=+/-`.

| computed block | symmetry-related partner | dimension |
|---|---|---:|
| `(0,0), C2=-` | self | 35,494 |
| `(0,0), C2=+` | self | 35,694 |
| `(0,1)` | `(0,3)` | 71,120 |
| `(0,2), C2=-` | self | 35,694 |
| `(0,2), C2=+` | self | 35,494 |
| `(1,0)` | `(2,0)` | 71,178 |
| `(1,1)` | `(2,3)` | 71,112 |
| `(1,2)` | `(2,2)` | 71,178 |
| `(1,3)` | `(2,1)` | 71,112 |

Star-weighting these blocks reconstructs all 853,776 states. The production
calculation retains six eigenvalues in each of the nine displayed blocks.

## Calculation and numerical checks

The complete scan covers `U/t=6:0.5:14`. Guarded block LOBPCG is used because a
direct benchmark found that multilevel Lanczos could return a ghost Ritz value.
All 153 coarse records passed explicit residual and orthogonality checks. The
largest explicit residual is `8.59e-9` and the largest orthogonality error is
`3.00e-15`.

The two ground-branch crossings were subsequently sampled every `0.05` in the
two relevant windows, using only `(0,0,C2=+)` and `(1,0)`. All 36 dense records
passed; their largest explicit residual is `1.49e-10`.

Independent validation includes exact geometry/group tests, reconstruction of a
complete six-site spectrum from momentum blocks, explicit fermionic-matrix
comparisons, all five `k <-> -k` partner checks, and reconstruction of the
translation-only spectra from both `C2` parities.

## Observed crossings

### Ground branch on `Ky=0`

Two crossings occur between different symmetry blocks,
`(0,0,C2=+)` and `(1,0)`. Because the states have different momentum quantum
numbers, these are symmetry-allowed exact crossings on this finite cluster.

| sampled crossing bracket | lower block to the left | lower block to the right |
|---|---|---|
| `6.60 < U/t < 6.65` | `(0,0,C2=+)` | `(1,0)` |
| `9.30 < U/t < 9.35` | `(1,0)` | `(0,0,C2=+)` |

Linear interpolation of the energy difference gives `U/t ~ 6.6281` and
`U/t ~ 9.3124`. These numbers are plot guides only; the actual ED resolution is
the `0.05` bracket.

### Other `Ky=0` crossings

The lowest `(0,0,C2=-)` and `(1,0)` branches also exchange order in the coarse
brackets `6.0--6.5` and `13.0--13.5`. The first is especially low lying (about
`0.007 t` above the absolute ground near the crossing) but is not itself an
absolute-ground crossing. Further crossings occur among higher retained levels;
they are visible in the dedicated `Ky=0` spectrum plot.

### `Ky=pi` crossings

There is no absolute-ground crossing on the `Ky=pi` line. There is, however, a
notably low crossing between the lowest `(0,2,C2=-)` and `(1,2)` branches in the
coarse bracket `9.5--10.0`. Both lie only about `0.01--0.02 t` above the absolute
ground at the bracket endpoints. Crossings among higher retained `Ky=pi` levels
also occur; the clearest coarse brackets include `8.0--8.5`, `8.5--9.0`,
`9.0--9.5`, `9.5--10.0`, `10.0--10.5`, `11.0--11.5`, `12.5--13.0`, and
`13.5--14.0`, depending on the pair of branches.

### Other low-energy rearrangements

The full set of block-lowest branches contains additional crossings involving
`(0,0,C2=-)`, `(0,1)`, `(0,2,C2=-)`, `(1,1)`, `(1,2)`, and `(1,3)`. Several
cluster near `U/t ~ 8.5--10`, where many distinct irreps lie close together.
This is a useful finite-size spectral observation, but the coarse scan alone does
not justify selecting every visual intersection for refinement.

Crossings between different listed irreps can be identified from energy order.
Within one irrep, sorted eigenvalues alone do not establish state identity across
`U/t`; overlaps or targeted observables would be needed to distinguish an avoided
crossing from a rapid rearrangement.

## Plots and machine-readable results

The generated plots are:

- `plots/N12/complete_inequivalent_spectrum.{pdf,png}`: six levels in every block;
- `plots/N12/all_sector_lowest_branches.{pdf,png}`: the lowest branch of every block;
- `plots/N12/Ky_0_spectrum.{pdf,png}`: six levels in the `Ky=0` blocks;
- `plots/N12/Ky_pi_spectrum.{pdf,png}`: six levels in the `Ky=pi` blocks;
- `plots/N12/Ky_0_pi_lowest_branches.{pdf,png}`: focused line comparison;
- `plots/N12/crossing_zoom.{pdf,png}`: dense energy differences at the two
  absolute-ground crossings.

The corresponding tables are `data/processed/N12/coarse_summary.csv` and
`dense_summary.csv`; JLD2 versions and plain-text analyses are stored beside them.
Raw records, processed tables, and generated plots are intentionally ignored by
Git because they are reproducible calculation products.

## Current conclusion and next diagnostic

The 12-site spectrum contains two resolved absolute-ground crossings on `Ky=0`,
one very low `Ky=pi` excited-state crossing, and several further crossings among
low block branches. The scientifically conservative next step is to choose which
of these finite-cluster candidates merit state-overlap or minimal-observable
diagnostics. A `4 x 4` calculation remains deferred.
