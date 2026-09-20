
# Prompt: Symmetry-Reduced ED of the \(3\times4\) Triangular-Lattice Hubbard Model with XDiag

I want to perform a small-cluster exact-diagonalization study of the half-filled triangular-lattice Hubbard model using **XDiag**, preferably through its Julia interface.

The purpose is to test whether the transition observed in iDMRG is associated with a low-energy level crossing or a strong rearrangement of the many-body spectrum.

The main calculation should be:

> Compute the lowest several eigenvalues in every **symmetry-inequivalent crystal-momentum sector** as a function of \(U/t\), with special attention to the \(K_y=0\) and \(K_y=\pi\) momentum lines.

Do **not** waste time diagonalizing momentum sectors related exactly by crystalline symmetry.

Do not include reduced-density-matrix, entanglement, or effective-spin-model calculations.

---

## 1. Hubbard model

Study

\[
H=
-t\sum_{\langle ij\rangle,\sigma}
\left(
c^\dagger_{i\sigma}c_{j\sigma}
+
c^\dagger_{j\sigma}c_{i\sigma}
\right)
+
U\sum_i n_{i\uparrow}n_{i\downarrow}.
\]

Set

\[
t=1,
\qquad
u=U/t.
\]

Use the half-filled sector

\[
N_\uparrow=N_\downarrow=6
\]

on \(N=12\) sites.

The raw fixed-particle-number Hilbert-space dimension is

\[
\binom{12}{6}^2=853\,776.
\]

---

## 2. Cluster geometry

Use triangular-lattice primitive vectors

\[
\mathbf a_1=(1,0),
\qquad
\mathbf a_2=
\left(
\frac12,\frac{\sqrt3}{2}
\right).
\]

Use the \(3\times4\) torus

\[
\boxed{
\mathbf T_1=3\mathbf a_1,
\qquad
\mathbf T_2=4\mathbf a_2.
}
\]

Equivalently,

```text
x = 0,1,2
y = 0,1,2,3
```

with periodicity

```text
(x+3,y) == (x,y)
(x,y+4) == (x,y)
```

and, for example,

```julia
site(x,y) = 1 + y + 4x
```

for zero-based \(x,y\).

Generate all geometry programmatically.

The three positive nearest-neighbor displacement vectors are

\[
\mathbf a_1,\qquad
\mathbf a_2,\qquad
\mathbf a_2-\mathbf a_1.
\]

Verify:

- 12 sites;
- coordination number 6;
- exactly 36 undirected nearest-neighbor bonds;
- no duplicated or missing bonds.

---

## 3. Translation symmetry and momentum labels

Define

\[
T_x:(x,y)\rightarrow(x+1,y),
\]

\[
T_y:(x,y)\rightarrow(x,y+1).
\]

They satisfy

\[
T_x^3=1,
\qquad
T_y^4=1.
\]

The allowed crystal momenta are labeled by

\[
(n_x,n_y),
\qquad
n_x=0,1,2,
\qquad
n_y=0,1,2,3,
\]

with

\[
T_x\rightarrow e^{i2\pi n_x/3},
\qquad
T_y\rightarrow e^{i2\pi n_y/4}.
\]

Thus

\[
K_y=\frac{2\pi n_y}{4}.
\]

Use `(nx, ny)` as the primary momentum label in the code.

---

## 4. Do not calculate all 12 momenta independently

For this \(3\times4\) torus, determine the surviving crystalline point group explicitly.

For the simple torus

\[
\mathbf T_1=(3,0),
\qquad
\mathbf T_2=(0,4),
\]

the nontrivial surviving point-group operation is the twofold rotation / inversion

\[
\boxed{
C_2:(x,y)\rightarrow(-x,-y).
}
\]

It sends

\[
(n_x,n_y)
\rightarrow
(-n_x\bmod3,\,-n_y\bmod4).
\]

Therefore momentum sectors related by

\[
\mathbf k\leftrightarrow-\mathbf k
\]

have identical spectra and only one representative needs to be diagonalized.

The 12 momenta reduce to the following **7 inequivalent momentum stars**:

| representative | symmetry-related partner |
|---|---|
| \((0,0)\) | self |
| \((0,1)\) | \((0,3)\) |
| \((0,2)\) | self |
| \((1,0)\) | \((2,0)\) |
| \((1,1)\) | \((2,3)\) |
| \((1,2)\) | \((2,2)\) |
| \((1,3)\) | \((2,1)\) |

Hence the production calculation should diagonalize only

```text
(0,0)
(0,1)
(0,2)
(1,0)
(1,1)
(1,2)
(1,3)
```

unless an explicit validation calculation requires the partner sector.

Before production runs, verify numerically that each paired sector has identical low-energy eigenvalues.

---

## 5. Special \(K_y=0\) and \(K_y=\pi\) lines

The main physics interest is in the two transverse-momentum lines

\[
K_y=0
\]

and

\[
K_y=\pi.
\]

### \(K_y=0\)

The three momenta are

```text
(0,0)
(1,0)
(2,0)
```

but

\[
(1,0)\sim(2,0)
\]

under \(C_2\).

Therefore only two inequivalent calculations are required:

\[
\boxed{
(0,0),\qquad(1,0).
}
\]

### \(K_y=\pi\)

The three momenta are

```text
(0,2)
(1,2)
(2,2)
```

but

\[
(1,2)\sim(2,2).
\]

Therefore only two inequivalent calculations are required:

\[
\boxed{
(0,2),\qquad(1,2).
}
\]

These four sectors should receive particular attention in the level-spectroscopy analysis.

The iDMRG ring-momentum labels \(K_y=0,\pi\) should **not** be assumed to be identical to the global ED momentum quantum number. The goal here is simply to inspect whether the finite Hubbard spectrum develops competing low-energy branches on these transverse-momentum lines.

---

## 6. Little-group crystalline symmetry

At generic momentum, \(C_2\) sends \(\mathbf k\) to a distinct momentum \(-\mathbf k\). It therefore relates two sectors but does not further split a single fixed-\(\mathbf k\) block.

The two momenta invariant under \(C_2\) are

\[
\boxed{
(0,0),\qquad(0,2).
}
\]

For these momenta,

\[
-\mathbf k=\mathbf k
\]

modulo a reciprocal lattice vector, so \(C_2\) belongs to the little group.

Therefore resolve these two sectors further into

\[
\boxed{C_2=+1,\qquad C_2=-1}
\]

whenever supported conveniently by XDiag.

Thus the most important symmetry labels are:

```text
(0,0,C2=+)
(0,0,C2=-)

(0,2,C2=+)
(0,2,C2=-)

(1,0)
(1,1)
(1,2)
(1,3)
(0,1)
```

Do not artificially combine \(+\mathbf k\) and \(-\mathbf k\) into a larger block if simply calculating one representative is cheaper.

---

## 7. XDiag symmetry implementation

Use XDiag `Permutation`, `PermutationGroup`, and `Representation` objects.

First construct the translation group

\[
Z_3\times Z_4.
\]

For each of the 7 inequivalent momenta, construct the corresponding one-dimensional translation representation.

Use

```julia
Electron(N, nup, ndn, representation)
```

with

```text
N   = 12
nup = 6
ndn = 6
```

for the symmetry-resolved Hubbard blocks.

For the self-conjugate momenta `(0,0)` and `(0,2)`, construct the appropriate enlarged little-group representations including \(C_2\) parity.

For each block print and save:

```text
nx
ny
Ky
C2_parity   # when applicable
block_dimension
```

---

## 8. Symmetry validation

Before the production scan:

1. verify that \(T_x\), \(T_y\), and \(C_2\) preserve all Hubbard hopping bonds;
2. verify the symmetry group algebra;
3. verify numerically that
   \[
   [H,T_x]=[H,T_y]=[H,C_2]=0;
   \]
4. verify
   \[
   E(n_x,n_y)=E(-n_x,-n_y);
   \]
5. explicitly compare the following pairs at several \(U/t\):

```text
(0,1) <-> (0,3)
(1,0) <-> (2,0)
(1,1) <-> (2,3)
(1,2) <-> (2,2)
(1,3) <-> (2,1)
```

6. once equivalence is confirmed, do not calculate both members during production scans;
7. check that the \(C_2=\pm1\) decomposition at `(0,0)` and `(0,2)` reproduces the translation-only spectrum.

---

## 9. Hamiltonian construction

Write

\[
H=H_t+UH_U
\]

with

\[
H_t=
-\sum_{\langle ij\rangle,\sigma}
\left(
c^\dagger_{i\sigma}c_{j\sigma}
+\mathrm{h.c.}
\right),
\]

\[
H_U=
\sum_i n_{i\uparrow}n_{i\downarrow}.
\]

Construct the Hamiltonian with XDiag `OpSum`.

Keep the geometry, symmetry construction, Hamiltonian, and analysis code separate.

Use Lanczos for the low-energy spectrum.

---

## 10. \(U/t\) scan

For every \(U/t\), calculate several lowest eigenvalues in the **7 inequivalent momentum stars only**.

For `(0,0)` and `(0,2)`, retain the additional \(C_2\) label when available.

Save

```text
U_over_t
energy
energy_minus_ground
nx
ny
Ky
C2_parity
level_index
block_dimension
Lanczos_residual
```

Use a coarse \(U/t\) scan first.

Once candidate crossings are found, perform a fine scan around them.

---

## 11. Main spectroscopy questions

Determine:

1. which momentum sector contains the absolute ground state as \(U/t\) changes;
2. whether the ground-state momentum changes;
3. whether two low-energy branches cross;
4. whether the crossing levels have different momentum quantum numbers;
5. whether levels at `(0,0)` or `(0,2)` have different \(C_2\) parity;
6. whether levels with identical symmetry quantum numbers show an avoided crossing;
7. whether the most important rearrangement occurs on the \(K_y=0\) and/or \(K_y=\pi\) lines;
8. whether its location is close to the transition observed in iDMRG.

---

## 12. Required plots

### Complete inequivalent spectrum

Plot the lowest several levels from the 7 inequivalent momentum stars versus \(U/t\).

### \(K_y=0\)

Plot

\[
(0,0),\qquad(1,0),
\]

including \(C_2\) parity for `(0,0)`.

### \(K_y=\pi\)

Plot

\[
(0,2),\qquad(1,2),
\]

including \(C_2\) parity for `(0,2)`.

### Combined \(K_y=0,\pi\)

Overlay the lowest branches from

```text
(0,0)
(1,0)
(0,2)
(1,2)
```

to expose possible crossings.

### Crossing zoom

After identifying a candidate crossing, perform a dense scan and generate a dedicated zoomed plot.

---

## 13. Observables only after identifying the crossing

Do not calculate many observables initially.

First understand the momentum-resolved spectrum.

For the few states directly involved in the crossing, then calculate:

### Spin structure factor

\[
S(\mathbf q)
=
\frac1N
\sum_{ij}
e^{i\mathbf q\cdot(\mathbf r_i-\mathbf r_j)}
\langle
\mathbf S_i\cdot\mathbf S_j
\rangle.
\]

### Scalar-chirality correlations

\[
\chi_{ijk}
=
\mathbf S_i\cdot
(\mathbf S_j\times\mathbf S_k).
\]

Calculate chirality-chirality correlations or a uniform chirality structure factor.

### Double occupancy

\[
D=
\frac1N
\sum_i
\langle n_{i\uparrow}n_{i\downarrow}\rangle.
\]

These are secondary diagnostics.

No reduced-density-matrix or entanglement calculation is required.

---

## 14. Development order

### Stage 1
Construct and validate the \(3\times4\) triangular torus.

### Stage 2
Implement the Hubbard Hamiltonian.

### Stage 3
Implement \(T_x,T_y\) and all momentum representations.

### Stage 4
Implement \(C_2\) and verify the momentum-star equivalences.

### Stage 5
Reduce the 12 momenta to the 7 inequivalent representatives.

### Stage 6
Implement \(C_2=\pm1\) little-group resolution at `(0,0)` and `(0,2)`.

### Stage 7
Validate spectra against the unreduced calculation at selected parameters.

### Stage 8
Perform the full \(U/t\) scan using only the inequivalent sectors.

### Stage 9
Analyze the \(K_y=0\) and \(K_y=\pi\) branches.

### Stage 10
Refine around candidate level crossings.

### Stage 11
Compute minimal physical observables for the crossing states.

The primary deliverable is a **symmetry-efficient momentum-resolved Hubbard spectrum versus \(U/t\)**, with special emphasis on the \(K_y=0\) and \(K_y=\pi\) momentum lines.
