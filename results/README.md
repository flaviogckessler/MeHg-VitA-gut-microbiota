# Results

- `alpha_diversity/` contains batch-specific alpha-diversity estimates and their source-associated metadata. Figure workflows join these estimates to the canonical manifests in `data/metadata/`.
- `taxon_relative_abundance/` contains batch-specific processed taxon relative-abundance matrices used by Figure 3 panels H-K and Figure 4 panels A-D.
- [`functional_prediction/`](functional_prediction/) contains validated PICRUSt2 EC, KO, MetaCyc, and weighted-NSTI outputs, separated by sequencing batch and renamed by content. Its `KO2` files remain explicitly marked as legacy because that transformation is undocumented. Each batch also contains generated `kegg_pathway_ko_z_scores/` sources for the Figure 7, Figure 8, and Figure S8-S16 heatmaps; the portable supplementary renderers and cohort mapping are documented in [`supplementary_figures/README.md`](../supplementary_figures/README.md).
- `legacy_summary_tables/` contains retained summary artifacts, including the cross-batch Z-score table, the Figure 4 `Signif_*.csv` validation checkpoints under [`legacy_summary_tables/figure_4/`](legacy_summary_tables/figure_4/), and the unchanged KO heatmap checkpoints under [`legacy_summary_tables/figure_7_8_s8_s16/`](legacy_summary_tables/figure_7_8_s8_s16/).

The `kegg_ortholog_abundance_legacy_ko2.tsv` filenames deliberately retain the `KO2` distinction because the precise transformation represented by `KO2` is not documented in the source files.

## Raw-sequence availability

The raw reads underlying these results are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample accessions span `SAMN61406324`–`SAMN61406377` for Dams_F0/F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365/F2-PND365. They are currently on hold and will become public upon publication of the manuscript.
