# Numerical data format

Numerical output is JLD2. Files contain only strings, numbers, arrays, and plain
string-keyed dictionaries; XDiag objects and project-specific structs are never
persisted.

## Raw records

Each raw file stores one `(symmetry block, U/t)` calculation under the JLD2 key
`result`. The schema version is 1. Main dictionaries are:

- `metadata`: scan name, configuration path, and UTC creation time;
- `cluster`: dimensions, particle numbers, and bond count;
- `parameters`: `U_over_t`, `t`, and `U`;
- `sector`: `(nx,ny)`, `Ky`, `C2_parity`, momentum star, and block dimension;
- `solver`: algorithm, convergence diagnostics, timings, seed, and CSR nonzeros;
- `spectrum`: energies in units of `t`, block-relative gaps, and level indices;
- `observables`: empty during the initial spectrum scan;
- `provenance`: Julia/XDiag versions, git revision, and Slurm identifiers.

`C2_parity=0` means that inversion is not in the momentum little group. It is not a
physical parity eigenvalue.

Raw records are written atomically through a temporary file. Existing results are
not overwritten by default.

## Processed summaries

The scan summary is stored under the JLD2 key `summary`. Its energy array is indexed
as

```text
energy[U_index, sector_index, level_index]
```

The summary computes `energy_minus_ground` relative to the absolute ground state
across all production blocks at each `U/t`; raw records only know their own
block-relative gaps. The summary builder rejects missing, unexpected, unconverged,
or metadata-inconsistent records.

The summary step also writes a long-form CSV beside the JLD2 file. It contains one
row per `(U/t, sector, level)` and is the portable plotting and table interface.

## Observable records and summaries

The coarse observable rerun uses the same raw schema, with `observables` containing
six-level arrays for double occupancy, translation-averaged spin correlations,
connected charge correlations, and the uniform chirality on each triangle
sublattice. Every rerun record stores its maximum energy difference from the
trusted coarse spectrum.

Ground-state chirality files additionally store the `2 x 2 x 3 x 4` complex Gram
matrix

```text
chirality_sublattice_gram[a, b, q_nx, q_ny]
```

for the up/down triangle-sublattice Fourier operators. Keeping this matrix is
required because the two triangle centers have different intra-cell offsets; a
single periodic `3 x 4` scalar array would discard their form factor.

The observable validator writes:

- `observables_summary_states.csv`: 918 rows, one per retained state;
- `observables_summary_ground.csv`: one ground-state row per `U/t`;
- `observables_summary_chirality_grid.csv`: first-Brillouin-zone chirality values
  at all 12 allowed momenta;
- `observables_summary_chirality_sublattice_gram.csv`: the sublattice Gram data
  used for the two representative heatmaps;
- `observables_summary_heatmap_correlations.csv`: spin and charge displacement
  correlations for those representative states.

The associated JLD2 summary records global residual, orthogonality, reference
energy, positivity, and uniform-chirality validation diagnostics.
