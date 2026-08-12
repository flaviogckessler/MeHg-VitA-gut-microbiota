# Alpha-diversity results

These batch-specific tables preserve the alpha-diversity estimates and associated source metadata used by the figure workflows.

- `batch01_iontorrent_dams-f0_f1-pnd21_alpha_diversity_with_metadata.csv` originates from `Alpha_filtered_Metadata_F0F1.csv`.
- `batch02_illumina-miseq_f1-f2_pnd365_alpha_diversity_with_metadata.csv` originates from `alpha_diversity.csv`.

Analyses should use the alpha-diversity values from these files and join sample attributes from the canonical manifests in `data/metadata/`. This excludes undocumented or non-primary samples without modifying the preserved result tables.

The tables are canonical numerical inputs for main Figure 3 panels D-G
(Shannon) and supplementary Figure S1 panels A-H (Observed and InvSimpson).
The Figure S1 workflow records the unmatched batch 02 row `t3` but excludes it
from analysis because no canonical metadata row was supplied.
