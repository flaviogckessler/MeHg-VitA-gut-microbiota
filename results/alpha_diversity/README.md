# Alpha-diversity results

These batch-specific tables preserve the alpha-diversity estimates and associated source metadata used by the figure workflows.

- `batch01_iontorrent_dams-f0_f1-pnd21_alpha_diversity_with_metadata.csv` originates from `Alpha_filtered_Metadata_F0F1.csv`.
- `batch02_illumina-miseq_f1-f2_pnd365_alpha_diversity_with_metadata.csv` originates from `alpha_diversity.csv`.

Analyses should use the alpha-diversity values from these files and join sample attributes from the canonical manifests in `data/metadata/`. This excludes undocumented or non-primary samples without modifying the preserved result tables.

The tables are canonical numerical inputs for main Figure 3 panels D-G
(Shannon) and supplementary Figure S1 panels A-H (Observed and InvSimpson).
The Figure S1 workflow records the unmatched batch 02 row `t3` but excludes it
from analysis because no canonical metadata row was supplied.

These R-derived microbiome summaries retain provenance to R and phyloseq [@r-core; @phyloseq].

## References

- <a id="ref-phyloseq"></a>McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE, 8*(4), e61217. [https://doi.org/10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
