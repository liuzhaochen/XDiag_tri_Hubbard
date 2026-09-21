# Symmetry-sector inventory

The cluster translations satisfy `Tx^3=1` and `Ty^4=1`. Momentum labels are
`nx=0,1,2` and `ny=0,1,2,3`, with `Ky=2*pi*ny/4`.

The point group is expected to reduce mechanically from triangular `D6` to
`{E,C2}`, where `C2: (x,y) -> (-x,-y)`. The one-dimensional irreps `A` and `B`
have `C2` characters `+1` and `-1`. The resulting production stars are:

| representative | partner | little-group resolution |
|---|---|---|
| `Gamma=(0,0)` | self | `A`, `B` |
| `(0,1)` | `(0,3)` | none |
| `M2=(0,2)` | self | `A`, `B` |
| `(1,0)` | `(2,0)` | none |
| `(1,1)` | `(2,3)` | none |
| `(1,2)` | `(2,2)` | none |
| `(1,3)` | `(2,1)` | none |

Validated with Julia 1.13.0 and XDiag 0.5.1:

| production block | star size | dimension |
|---|---:|---:|
| `Gamma(B)` | 1 | 35,494 |
| `Gamma(A)` | 1 | 35,694 |
| `(0,1)` | 2 | 71,120 |
| `M2(B)` | 1 | 35,694 |
| `M2(A)` | 1 | 35,494 |
| `(1,0)` | 2 | 71,178 |
| `(1,1)` | 2 | 71,112 |
| `(1,2)` | 2 | 71,178 |
| `(1,3)` | 2 | 71,112 |

The two irrep dimensions sum to 71,188 at both self-conjugate momenta. Weighting
translation-only representative dimensions by star size reconstructs all 853,776
fixed-particle-number states.
