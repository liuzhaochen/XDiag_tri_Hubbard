# Spectrum plots

After the validated coarse summary and long-form CSV have been created, generate the
four required overview figures with:

```bash
python3 scripts/plot_spectrum.py \
  data/processed/N12/coarse_summary.csv plots/N12
```

The script produces PDF and PNG versions of the complete inequivalent spectrum, the
`Ky=0` spectrum, the `Ky=pi` spectrum, and the combined lowest branches. A dedicated
crossing zoom will be added only after the coarse data locate a candidate crossing.
