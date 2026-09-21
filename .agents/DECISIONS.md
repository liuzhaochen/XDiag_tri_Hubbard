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
- `(0,0)` and `(0,2)` are additionally resolved by `C2` parity.
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
  `(0,0,C2=+)` and `(1,0)` blocks. Other low-lying crossings remain explicit
  candidates for discussion before further targeted work.
- The reader-facing scientific narrative is maintained in
  `docs/physical_report.md`; `.agents` remains the compact operational ledger.
- Treat the intermediate ground as a twofold `(1,0)/(2,0)` momentum-star pair.
  The label has `Tx=exp(+/- i 2 pi/3)` and `Ty=1` and is not identified with the
  conventional two-dimensional triangular-lattice `K` point.
