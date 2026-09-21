# Spectrum plots

After the validated coarse summary and long-form CSV have been created, generate the
overview figures with:

```bash
python3 scripts/plot_spectrum.py \
  data/processed/N12/coarse_summary.csv plots/N12
```

The script produces PDF and PNG versions of the complete inequivalent spectrum, the
lowest branch from every block, the `Ky=0` spectrum, the `Ky=pi` spectrum, and the
combined `Ky=0,pi` lowest branches. Generate the refined two-crossing zoom with:

```bash
python3 scripts/plot_crossing_zoom.py \
  data/processed/N12/dense_summary.csv plots/N12/crossing_zoom
```

Interpretation and a complete plot index are in `docs/physical_report.md`.

Generate the cluster geometry and allowed-momentum schematic with:

```bash
python3 scripts/plot_cluster_schematic.py plots/N12
```

After the observable validator creates the five processed observable tables,
generate the ground-state, low-state, and structure-factor figures with:

```bash
python3 scripts/plot_observables.py \
  data/processed/N12/observables_summary_states.csv \
  data/processed/N12/observables_summary_ground.csv \
  data/processed/N12/observables_summary_heatmap_correlations.csv \
  data/processed/N12/observables_summary_chirality_grid.csv \
  data/processed/N12/observables_summary_chirality_sublattice_gram.csv \
  plots/N12
```

The heatmap colors are lower-resolution finite Fourier interpolations over an
extended momentum region. White dots show allowed-momentum replicas and the white
dashed hexagon marks the first Brillouin-zone boundary.
