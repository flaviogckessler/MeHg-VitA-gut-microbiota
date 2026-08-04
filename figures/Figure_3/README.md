# Figure 3 panels A-C

The scripts use canonical repository data directly. Figure-specific copies are not required.

## Canonical inputs

- `data/processed/batch01_iontorrent_dams-f0_f1-pnd21/asv_counts_with_taxonomy.tsv`
- `data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv`
- `data/processed/batch02_illumina-miseq_f1-f2_pnd365/asv_relative_abundance_with_taxonomy.tsv`
- `data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv`

## Run from the repository root

```r
source("figures/Figure_3/validate_figure_3_abc.R")
source("figures/Figure_3/generate_figure_3_abc.R")
```

Generated panels are written as lightweight PNG files in `figures/Figure_3/output`.
