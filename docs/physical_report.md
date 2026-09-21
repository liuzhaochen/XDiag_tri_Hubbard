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
representative of each star. Only the self-conjugate momenta `Gamma=(0,0)` and
`M2=(0,2)` have the little co-group `C2`. We label its one-dimensional irreps
`A` and `B`, with `C2` characters `+1` and `-1`, respectively. These are irreps
of the surviving cluster little group, not irreps of the full triangular `D6`.

| computed block | symmetry-related partner | dimension |
|---|---|---:|
| `Gamma(B)` | self | 35,494 |
| `Gamma(A)` | self | 35,694 |
| `(0,1)` | `(0,3)` | 71,120 |
| `M2(B)` | self | 35,694 |
| `M2(A)` | self | 35,494 |
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
two relevant windows, using only `Gamma(A)` and `(1,0)`. All 36 dense records
passed; their largest explicit residual is `1.49e-10`.

Independent validation includes exact geometry/group tests, reconstruction of a
complete six-site spectrum from momentum blocks, explicit fermionic-matrix
comparisons, all five `k <-> -k` partner checks, and reconstruction of the
translation-only spectra from both `A/B` little-group irreps.

The observable rerun recomputed all six retained states in all 153 coarse-grid
blocks. Every spectrum agrees with the original production energy to within
`3.02e-14 t`; the largest explicit residual remains `8.59e-9`. A separate
17-task calculation measured the ground-state chirality structure factor. In
total, 170 raw observable records and 918 state rows pass the validator.

## Observable definitions and momentum convention

For each retained state we measure the double occupancy per site,

\[
D=\frac{1}{N}\sum_i\langle n_{i\uparrow}n_{i\downarrow}\rangle,
\]

and translation-averaged spin and connected charge structure factors,

\[
S(\mathbf q)=\sum_{\mathbf r}e^{i\mathbf q\cdot\mathbf r}
\langle\mathbf S_0\cdot\mathbf S_{\mathbf r}\rangle,
\qquad
N(\mathbf q)=\sum_{\mathbf r}e^{i\mathbf q\cdot\mathbf r}
\langle(n_0-1)(n_{\mathbf r}-1)\rangle.
\]

The conventional triangular-lattice corner `K=(-1/3,1/3)` is not an allowed
momentum of this torus. Plots use the closest allowed point
`Kcl=(-1/3,1/4)`, whose translation-sector label is `(nx,ny)=(2,1)`. Of the
three `M` representatives, `M2=(0,-1/2)` is exactly allowed; the other two are
represented by `M1cl=(1/3,0)` and `M3cl=(1/3,1/2)`. Their sector labels are
`(1,0)`, `(0,2)`, and `(1,2)`, respectively. The cluster schematic explicitly
shows all 12 allowed momenta and these choices.

Scalar chirality is

\[
\chi_{ijk}=\mathbf S_i\cdot(\mathbf S_j\times\mathbf S_k)
=i[\mathbf S_i\cdot\mathbf S_j,\mathbf S_j\cdot\mathbf S_k],
\]

with all 24 elementary up/down triangles oriented counterclockwise. The
connected structure factor is normalized per triangle:

\[
C_\chi(\mathbf q)=\frac{1}{N_\triangle}
\left\langle
\left(\sum_\triangle e^{i\mathbf q\cdot\mathbf r_\triangle}\chi_\triangle\right)^\dagger
\left(\sum_\triangle e^{i\mathbf q\cdot\mathbf r_\triangle}\chi_\triangle\right)
\right\rangle
-N_\triangle|\bar\chi|^2\delta_{\mathbf q,0},
\quad
\bar\chi=\frac{1}{N_\triangle}\sum_\triangle\langle\chi_\triangle\rangle,
\]

where `Ntriangle=24`. The subtraction affects only `Gamma`. The implementation
retains the full `2 x 2` Gram matrix of the up/down-triangle Fourier operators,
including their distinct intra-cell centers. This is necessary to fold momenta
and make the continuous momentum-space visualization without dropping the
two-sublattice form factor. Each allowed-momentum Gram matrix is checked to be
Hermitian and positive semidefinite.

For the finite-momentum ground region, the degenerate ground space is the
time-reversal pair `{|k>,T|k>}` in `(1,0)` and `(2,0)`. Uniform chirality cannot
mix the two states because it preserves momentum, while time reversal makes
their diagonal matrix elements opposite. Thus its projected eigenvalues are
exactly `+/- <k|chibar|k>`, which are evaluated from one representative. At
`Gamma` the sampled ground is nondegenerate, so the projection is one
dimensional.

## Observed crossings

### Ground branch on `Ky=0`

Two crossings occur between different symmetry blocks,
`Gamma(A)` and `(1,0)`. Because the states have different momentum quantum
numbers, these are symmetry-allowed exact crossings on this finite cluster.

| sampled crossing bracket | lower block to the left | lower block to the right |
|---|---|---|
| `6.60 < U/t < 6.65` | `Gamma(A)` | `(1,0)` |
| `9.30 < U/t < 9.35` | `(1,0)` | `Gamma(A)` |

Linear interpolation of the energy difference gives `U/t ~ 6.6281` and
`U/t ~ 9.3124`. These numbers are plot guides only; the actual ED resolution is
the `0.05` bracket.

An additional audit diagonalized all 12 translation momenta without `C2`
resolution, using a new random seed and explicit residual checks. It independently
reproduced both order reversals:

| `U/t` | `E(0,0)` | `E(1,0)=E(2,0)` | `E(1,0)-E(0,0)` |
|---:|---:|---:|---:|
| 6.60 | -8.595121557329 | -8.594593337081 | +0.000528220248 |
| 6.65 | -8.527904893136 | -8.528316662264 | -0.000411769128 |
| 8.00 | -6.957821307096 | -6.991271041267 | -0.033449734171 |
| 9.30 | -5.935517088312 | -5.936402816482 | -0.000885728169 |
| 9.35 | -5.905634931836 | -5.902940274534 | +0.002694657302 |

The `(1,0)` and `(2,0)` partner energies agree within `1.6e-14 t`; the largest
explicit residual in this audit is `5.51e-9`. At `U/t=8`, a separate all-momentum
`(Nup,Ndn)=(7,5)` calculation has minimum energy `-6.913712092311 t`, which is
`0.077558948956 t` above the `(6,6)` ground. By spin SU(2), the finite-momentum
ground is therefore consistent with total spin `S=0`, rather than an unnoticed
higher-spin multiplet.

The finite-momentum ground is a degenerate `k`/`-k` pair, not a unique complex
state. Here `(1,0)` means translation eigenvalues
`Tx=exp(i 2 pi/3)` and `Ty=1`; it should not be conflated with the conventional
two-dimensional triangular-lattice `K` point. Such a momentum-changing ground
level is allowed on this non-bipartite, low-symmetry finite torus.

### Other `Ky=0` crossings

The lowest `Gamma(B)` and `(1,0)` branches also exchange order in the coarse
brackets `6.0--6.5` and `13.0--13.5`. The first is especially low lying (about
`0.007 t` above the absolute ground near the crossing) but is not itself an
absolute-ground crossing. Further crossings occur among higher retained levels;
they are visible in the dedicated `Ky=0` spectrum plot.

### `Ky=pi` crossings

There is no absolute-ground crossing on the `Ky=pi` line. There is, however, a
notably low crossing between the lowest `M2(B)` and `(1,2)` branches in the
coarse bracket `9.5--10.0`. Both lie only about `0.01--0.02 t` above the absolute
ground at the bracket endpoints. Crossings among higher retained `Ky=pi` levels
also occur; the clearest coarse brackets include `8.0--8.5`, `8.5--9.0`,
`9.0--9.5`, `9.5--10.0`, `10.0--10.5`, `11.0--11.5`, `12.5--13.0`, and
`13.5--14.0`, depending on the pair of branches.

### Other low-energy rearrangements

The full set of block-lowest branches contains additional crossings involving
`Gamma(B)`, `(0,1)`, `M2(B)`, `(1,1)`, `(1,2)`, and `(1,3)`. Several
cluster near `U/t ~ 8.5--10`, where many distinct irreps lie close together.
This is a useful finite-size spectral observation, but the coarse scan alone does
not justify selecting every visual intersection for refinement.

Crossings between different listed irreps can be identified from energy order.
Within one irrep, sorted eigenvalues alone do not establish state identity across
`U/t`; overlaps or targeted observables would be needed to distinguish an avoided
crossing from a rapid rearrangement.

## Observable results

### Ground branch

The double occupancy decreases from `0.124682` at `U/t=6` to `0.021494` at
`U/t=14`. The first momentum-changing ground crossing produces only a modest
observable change. The second produces a much clearer branch switch: the
`Gamma(A)` branch has stronger `Kcl` spin correlations and weaker `M2`
correlations than the `(1,0)` branch.

An earlier figure used `mean Mcl`, the arithmetic mean of the three values at
`M1cl`, exact `M2`, and `M3cl`, with a shaded minimum-to-maximum range. Because
the `3 x 4` torus does not preserve the symmetry relating the three continuum
`M` points, that average hid meaningful cluster anisotropy. The revised figure
shows all three curves separately.

| `U/t` | ground block | `D` | `S(Kcl)` | `S(M2)` | `N(Kcl)` |
|---:|---|---:|---:|---:|---:|
| 6.0 | `Gamma(A)` | 0.124682 | 0.811082 | 0.756704 | 0.300860 |
| 8.0 | `(1,0)` | 0.080198 | 1.017817 | 1.067291 | 0.207080 |
| 9.0 | `(1,0)` | 0.061089 | 1.123415 | 1.199851 | 0.161851 |
| 9.5 | `Gamma(A)` | 0.046730 | 1.548072 | 1.031221 | 0.130959 |
| 14.0 | `Gamma(A)` | 0.021494 | 2.016901 | 0.778196 | 0.065011 |

At `U/t=8`, `S(Kcl)` and `S(M2)` remain comparable. At `U/t=14`, the much
larger `S(Kcl)` is consistent with increasingly strong finite-cluster
120-degree-like spin correlations. This statement concerns the structure factor
on the 12-site torus, not an order-parameter extrapolation.

The charge structure factors decrease smoothly as charge fluctuations are
suppressed. Both `S(Gamma)` and connected `N(Gamma)` vanish to numerical
precision in the singlet, fixed-particle-number ground states.

### Low-lying states and the `Ky=pi` crossing

At `U/t=8`, the retained `(1,0)` levels have gaps
`0, 0.57244, 0.59012, 0.66484, 0.78806, 0.88425 t`. Thus there is only one
state in the representative `(1,0)` block below `0.5 t`, but it has the exact
time-reversal partner in `(2,0)`. Globally, several other irreps lie below
`0.3 t`, as shown in the low-state observable plot.

The low `Ky=pi` crossing is also visible in observables. At `U/t=9.5`, the
lowest `(1,2)` state lies `0.013858 t` above the ground and just below the
`M2(B)` state at `0.014963 t`; at `U/t=10`, their gaps are respectively
`0.023205 t` and `0.010907 t`. The states have clearly different spin patterns:
at `U/t=10`, `(S(Kcl),S(M2))` is `(0.9103,0.8264)` for `(1,2)` and
`(1.4638,1.7243)` for `M2(B)`.

### Chirality diagnostics and the CSL question

No sampled ground state has a resolvable uniform scalar chirality. Projection is
only a useful statement in the degenerate `(1,0)/(2,0)` ground space, where the
two eigenvalues have magnitude at most `2.1e-13`. For an isolated `Gamma(A)`
singlet there is no nontrivial projected subspace: it is simply the ordinary
expectation value, which vanishes by symmetry. For this reason projected
`|chibar|` is no longer included as a ground-state plot curve. Across all 918
retained symmetry eigenstates, the largest absolute uniform chirality expectation
is `1.33e-11`, at the scale of numerical noise.

The connected chirality structure factor is nonzero because it measures the
variance and spatial correlations of chirality, not its expectation value.
Representative fixed-momentum values are:

| `U/t` | ground block | `Cchi(Gamma)` | `Cchi(Kcl)` | `Cchi(M2)` |
|---:|---|---:|---:|---:|
| 6.0 | `Gamma(A)` | 0.085970 | 0.055455 | 0.085365 |
| 8.0 | `(1,0)` | 0.078778 | 0.074439 | 0.080555 |
| 9.0 | `(1,0)` | 0.089739 | 0.078331 | 0.080416 |
| 9.5 | `Gamma(A)` | 0.264839 | 0.076022 | 0.166102 |
| 14.0 | `Gamma(A)` | 0.328275 | 0.065763 | 0.193415 |

The growth of `Cchi(Gamma)` above the second crossing does not imply a finite
chiral order parameter. In the isolated time-reversal-even `Gamma(A)` singlet,
`<chibar>=0`, while `Cchi(Gamma)` is proportional to
`<chibar^2>-<chibar>^2` and may be large. Increasing `U/t` suppresses double
occupancy and forms stronger local moments, which can enhance multi-spin
fluctuations even when time reversal is unbroken. Establishing chiral order would
require size scaling of the connected peak together with a nonzero projected
order parameter in a symmetry-breaking manifold. The former is unavailable on
one cluster and the latter is absent here. The displayed `U/t=8`, `(1,0)` map is
state-specific; its time-reversal partner is reflected under `q -> -q`.

These data do **not** justify identifying the `(1,0)` ground region as a chiral
spin liquid. A momentum doublet and finite chirality fluctuations are not, by
themselves, topological or chiral-order diagnostics; here the projected uniform
chirality is numerically zero. Establishing a CSL would require additional
evidence such as a robust quasi-degenerate topological manifold across cluster
geometries or twists, flux response/spectral flow, and consistent entanglement
or modular diagnostics. The present result should remain labeled a
finite-cluster momentum-changing ground branch.

## Plots and machine-readable results

The generated plots are:

- `plots/N12/complete_inequivalent_spectrum.{pdf,png}`: six levels in every block;
- `plots/N12/all_sector_lowest_branches.{pdf,png}`: the lowest branch of every block;
- `plots/N12/Ky_0_spectrum.{pdf,png}`: six levels in the `Ky=0` blocks;
- `plots/N12/Ky_pi_spectrum.{pdf,png}`: six levels in the `Ky=pi` blocks;
- `plots/N12/Ky_0_pi_lowest_branches.{pdf,png}`: focused line comparison;
- `plots/N12/crossing_zoom.{pdf,png}`: dense energy differences at the two
  absolute-ground crossings.
- `plots/N12/cluster_schematic.{pdf,png}`: real-space torus and all 12 allowed
  momentum points, including the selected `Kcl`;
- `plots/N12/ground_state_observables.{pdf,png}`: `D`, selected spin/charge
  structure factors, and chirality diagnostics over the 17-point grid;
- `plots/N12/low_lying_state_observables.{pdf,png}`: all retained states within
  `0.6 t` of the ground, colored by observables;
- `plots/N12/structure_factor_heatmaps.{pdf,png}`: lower-resolution extended-zone
  Fourier visualizations for `U/t=8`, `(1,0)` and `U/t=14`, `Gamma(A)`. White
  dots are allowed-momentum replicas and the white dashed hexagon is the first-BZ
  boundary. Between the dots the color field is a finite Fourier interpolation,
  not additional ED momentum resolution.

The corresponding spectrum tables are `coarse_summary.csv` and
`dense_summary.csv`. Observable products are `observables_summary_states.csv`,
`observables_summary_ground.csv`, `observables_summary_chirality_grid.csv`,
`observables_summary_chirality_sublattice_gram.csv`, and
`observables_summary_heatmap_correlations.csv` under `data/processed/N12/`.
JLD2 summaries are stored beside them. Raw records, processed tables, and
generated plots are intentionally ignored by Git because they are reproducible
calculation products.

## Current conclusion and next diagnostic

The 12-site spectrum contains two resolved absolute-ground crossings on `Ky=0`,
one very low `Ky=pi` excited-state crossing, and several further crossings among
low block branches. The complete coarse-grid observable rerun characterizes the
states on both sides without providing evidence for spontaneous scalar chirality
in the intermediate `(1,0)/(2,0)` ground doublet. The next discussion can decide
whether overlaps/twisted-boundary spectral flow are worthwhile and whether any
`4 x 4` calculation is scientifically and computationally justified; `4 x 4`
remains deferred for now.
