# Results

- `alpha_diversity/` contains batch-specific alpha-diversity estimates and their source-associated metadata. Figure workflows join these estimates to the canonical manifests in `data/metadata/`.
- `taxon_relative_abundance/` contains batch-specific processed taxon relative-abundance matrices used by Figure 3 panels H-K and Figure 4 panels A-D.
- [`functional_prediction/`](functional_prediction/) contains validated PICRUSt2 EC, KO, MetaCyc, and weighted-NSTI outputs, separated by sequencing batch and renamed by content. Its `KO2` files remain explicitly marked as legacy because that transformation is undocumented. Each batch also contains generated `kegg_pathway_ko_z_scores/` sources for the Figure 7, Figure 8, and Figure S8-S16 heatmaps; the portable supplementary renderers and cohort mapping are documented in [`supplementary_figures/README.md`](../supplementary_figures/README.md).
- `legacy_summary_tables/` contains retained summary artifacts, including the cross-batch Z-score table, the Figure 4 `Signif_*.csv` validation checkpoints under [`legacy_summary_tables/figure_4/`](legacy_summary_tables/figure_4/), and the unchanged KO heatmap checkpoints under [`legacy_summary_tables/figure_7_8_s8_s16/`](legacy_summary_tables/figure_7_8_s8_s16/).

The `kegg_ortholog_abundance_legacy_ko2.tsv` filenames deliberately retain the `KO2` distinction because the precise transformation represented by `KO2` is not documented in the source files.

## Raw-sequence availability

The raw reads underlying these results are publicly available in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample accessions span `SAMN61406324`–`SAMN61406377` for Dams_F0/F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365/F2-PND365.

Relevant methods and repositories are PICRUSt2, KEGG, R, and the Sequence Read Archive [@picrust2; @kegg; @r-core; @sra].

## References

- <a id="ref-picrust2"></a>Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., Huttenhower, C., & Langille, M. G. I. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology, 38*(6), 685–688. [https://doi.org/10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
- <a id="ref-kegg"></a>Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27–30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- <a id="ref-sra"></a>Leinonen, R., Sugawara, H., Shumway, M., & International Nucleotide Sequence Database Collaboration. (2011). The Sequence Read Archive. *Nucleic Acids Research, 39*(Database issue), D19–D21. [https://doi.org/10.1093/nar/gkq1019](https://doi.org/10.1093/nar/gkq1019)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
