# Symmetry-sector inventory

The cluster translations satisfy `Tx^3=1` and `Ty^4=1`. Momentum labels are
`nx=0,1,2` and `ny=0,1,2,3`, with `Ky=2*pi*ny/4`.

The point group is expected to reduce mechanically from triangular `D6` to
`{E,C2}`, where `C2: (x,y) -> (-x,-y)`. The resulting production stars are:

| representative | partner | little-group resolution |
|---|---|---|
| `(0,0)` | self | `C2=+`, `C2=-` |
| `(0,1)` | `(0,3)` | none |
| `(0,2)` | self | `C2=+`, `C2=-` |
| `(1,0)` | `(2,0)` | none |
| `(1,1)` | `(2,3)` | none |
| `(1,2)` | `(2,2)` | none |
| `(1,3)` | `(2,1)` | none |

Validated with Julia 1.13.0 and XDiag 0.5.1:

| production block | star size | dimension |
|---|---:|---:|
| `(0,0), C2=-` | 1 | 35,494 |
| `(0,0), C2=+` | 1 | 35,694 |
| `(0,1)` | 2 | 71,120 |
| `(0,2), C2=-` | 1 | 35,694 |
| `(0,2), C2=+` | 1 | 35,494 |
| `(1,0)` | 2 | 71,178 |
| `(1,1)` | 2 | 71,112 |
| `(1,2)` | 2 | 71,178 |
| `(1,3)` | 2 | 71,112 |

The two parity dimensions sum to 71,188 at both self-conjugate momenta. Weighting
translation-only representative dimensions by star size reconstructs all 853,776
fixed-particle-number states.
