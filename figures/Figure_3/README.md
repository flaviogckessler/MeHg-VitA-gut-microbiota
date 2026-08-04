# Figure 3

The scripts use canonical repository data directly. Figure-specific copies are not required.

## Canonical inputs for panels A-C

- `data/processed/batch01_iontorrent_dams-f0_f1-pnd21/asv_counts_with_taxonomy.tsv`
- `data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv`
- `data/processed/batch02_illumina-miseq_f1-f2_pnd365/asv_relative_abundance_with_taxonomy.tsv`
- `data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv`

## Canonical inputs for panels D-G

- `results/alpha_diversity/batch01_iontorrent_dams-f0_f1-pnd21_alpha_diversity_with_metadata.csv`
- `results/alpha_diversity/batch02_illumina-miseq_f1-f2_pnd365_alpha_diversity_with_metadata.csv`
- `data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv`
- `data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv`

## Canonical inputs for panels H-K

- `results/taxon_relative_abundance/batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv`
- `results/taxon_relative_abundance/batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv`
- `data/metadata/batch01_iontorrent_dams-f0_f1-pnd21/sample_manifest.csv`
- `data/metadata/batch02_illumina-miseq_f1-f2_pnd365/sample_manifest.csv`

## Run from the repository root

```r
source("figures/Figure_3/validate_figure_3_abc.R")
source("figures/Figure_3/generate_figure_3_abc.R")
source("figures/Figure_3/validate_figure_3_dg.R")
source("figures/Figure_3/generate_figure_3_dg.R")
source("figures/Figure_3/validate_figure_3_hk.R")
source("figures/Figure_3/generate_figure_3_hk.R")
```

Generated panels are written as lightweight PNG files in `figures/Figure_3/output`.
