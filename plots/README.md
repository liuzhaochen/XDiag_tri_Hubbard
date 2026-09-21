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
