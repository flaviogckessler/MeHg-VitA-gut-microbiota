# Taxon relative abundance

- `batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv` preserves the values from the original source file named `AbundanciasTaxons_Todos.csv`.
- `batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv` preserves the values from the original source file named `AbundanciasTaxons_1year.csv`.

Each table contains taxonomic descriptors in rows and sample-level relative abundances in columns. These files are the canonical repository locations used at runtime; figure workflows select the required taxa and join sample identifiers to the canonical manifests in [`data/metadata/`](../../data/metadata/).

The historical microbiome processing provenance uses R and phyloseq [@r-core; @phyloseq].

## References

- <a id="ref-phyloseq"></a>McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE, 8*(4), e61217. [https://doi.org/10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
