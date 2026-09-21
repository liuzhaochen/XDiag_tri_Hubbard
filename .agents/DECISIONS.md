# Decision log

## 2026-09-20

- Julia and XDiag are the implementation stack. Julia is never run on the login node.
- Tests and small calculations use Slurm's `test` partition, whose walltime limit is
  one hour.
- The initial and only active cluster is the `3 x 4` torus with 12 sites and
  `Nup=Ndn=6`. A `4 x 4` calculation is deferred until the 12-site study is complete.
- Site coordinates are triangular primitive-vector coefficients, with
  `site(x,y) = 1 + y + 4x` for zero-based `x=0:2`, `y=0:3`.
- The surviving point group must be derived from triangular-lattice `D6`, not inserted
  as a rectangular-lattice assumption.
- Production uses the seven momentum-star representatives specified in `prompt.md`.
- `Gamma=(0,0)` and `M2=(0,2)` are additionally resolved into the `A/B`
  irreps of their `C2` little group.
- The initial coarse grid is `U/t = 6:0.5:14`, with six physical levels requested
  from every block.
- Numerical results will use portable plain-container JLD2 records and atomic writes.
- Wavefunctions are temporary during the initial spectrum scan.
- Requested git commits use author `zhaochen liu <zhaochen.liu@uni-wuerzburg.de>`.
- Multilevel production uses CSR32 block LOBPCG with two guard vectors, explicit
  residual/orthogonality gates, and deterministic seed retries. On the largest
  `(1,0)` block at `U/t=8`, XDiag multilevel Lanczos reported convergence but
  returned a ghost Ritz value, residuals up to `1.11e-1`, and orthogonality error
  `9.996e-1`. LOBPCG returned six orthogonal states with residuals below `1e-10`.
- Single-level Lanczos remains useful as an independent ground-energy check; an
  unreliable multilevel result is never accepted merely because its criterion says
  `converged`.
- Production keeps both gates: XDiag must report `converged` and the explicit
  residual must be below `1e-8`. The sole coarse-array outlier `(1,3), U/t=14`
  reached `8.93e-9` at 500 iterations but was rejected because its criterion remained
  `maxiterations`; its targeted rerun uses a 1000-iteration cap.

## 2026-09-21

- Interpret the 12-site ED only as a finite-cluster crossing diagnostic. Neither
  crossing brackets nor linear interpolation are presented as thermodynamic phase
  boundaries or as a reproduction of iDMRG.
- The two absolute-ground crossings are refined only in the competing
  `Gamma(A)` and `(1,0)` blocks. Other low-lying crossings remain explicit
  candidates for discussion before further targeted work.
- The reader-facing scientific narrative is maintained in
  `docs/physical_report.md`; `.agents` remains the compact operational ledger.
- Treat the intermediate ground as a twofold `(1,0)/(2,0)` momentum-star pair.
  The label has `Tx=exp(+/- i 2 pi/3)` and `Ty=1` and is not identified with the
  conventional two-dimensional triangular-lattice `K` point.
- The observable rerun covers the existing 17-point coarse grid only; no new dense
  observable scan is needed because the crossing brackets are already resolved.
- Use `Kcl=(-1/3,1/4)`, sector label `(2,1)`, as the closest allowed momentum to
  the conventional triangular `K=(-1/3,1/3)` point.
- Scalar chirality uses all 24 counterclockwise elementary triangles. Its momentum
  structure retains the full up/down-triangle `2 x 2` Gram matrix so reciprocal
  folding preserves the intra-cell form factor.
- A continuous heatmap is a finite Fourier interpolation only. Every such plot
  overlays the 12 actually allowed momenta.
- Do not label the `(1,0)/(2,0)` region a CSL from this cluster: projected uniform
  chirality is numerically zero, and finite chirality fluctuations plus a momentum
  doublet are insufficient CSL diagnostics.
- Present self-conjugate sectors by professional little-group irreps:
  `Gamma(A/B)` and `M2(A/B)`, where `A/B` have `C2` characters `+1/-1`.
- Ground-state plots show the three cluster `M` representatives separately;
  do not average them into an ambiguous `mean Mcl` on this anisotropic torus.
- Do not plot projected uniform chirality for isolated singlets, and do not use
  `max_q Cchi(q)` as an order diagnostic. Show fixed-momentum connected chirality
  fluctuations at `Gamma`, `Kcl`, and `M2` instead.
- Momentum heatmaps use a lower-resolution extended-zone view with allowed-point
  replicas and a dashed first-BZ boundary; the heatmap mesh is rasterized in PDF.
- The validated `3 x 4` calculation is the completed active study. A possible
  `4 x 4` extension is deferred until the user requests it.
- The `4 x 4` torus restores the complete triangular `D6` point group. Its sector
  code must therefore implement the appropriate `D6` little-group irreps; the
  `3 x 4` `C2`-only labels must not be carried over.
- A translation-only `4 x 4` half-filled block is estimated to have dimension
  about `10,352,306`, about `540 million` CSR nonzeros, and about `10.9 GB` of
  CSR storage. Use a compute-node pilot before choosing production resources.
