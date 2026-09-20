# Data directory

- `raw/N12/<scan>/<sector>/u_<milli>.jld2`: one atomic calculation record per
  symmetry block and coupling.
- `processed/N12/`: validated scan summaries used for analysis and plotting.

Raw and processed JLD2 files are intentionally ignored by git. Their portable schema
is documented in [`../docs/data_format.md`](../docs/data_format.md). Compact reports,
tables, and figure source data may be tracked once the scan is complete.
