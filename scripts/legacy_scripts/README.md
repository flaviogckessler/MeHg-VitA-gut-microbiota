# Legacy scripts used for analysis and figure provenance

These historical scripts were used as references when portable analysis and figure workflows were created. Their contents are preserved unchanged for provenance; they may contain historical paths, interactive state assumptions, package-installation calls, overwritten exploratory commands, and obsolete code, so validated workflows do not source or execute them at runtime.

| Legacy reference script | Portable workflow derived from it |
|---|---|
| [`dada2_phyloseq_batch01_iontorrent.R`](dada2_phyloseq_batch01_iontorrent.R) | [`preprocess_16s.R`](../01_sequence_processing/preprocess_16s.R), batch 01 single-end Ion Torrent filtering, denoising, chimera removal, and processed-data exports. |
| [`dada2_batch02_illumina.R`](dada2_batch02_illumina.R) | [`preprocess_16s.R`](../01_sequence_processing/preprocess_16s.R), batch 02 paired-end Illumina filtering, denoising, merging, chimera removal, and processed-data exports. |
| [`phyloseq_exploratory.R`](phyloseq_exploratory.R) and [`taxa_summary.R`](taxa_summary.R) | Preserved as downstream exploratory provenance. They are not part of FASTQ preprocessing and are not sourced by the canonical workflow. |
| [`Abundancias_Barras_Mestrado.R`](Abundancias_Barras_Mestrado.R) | [`generate_figure_3_abc.R`](../../figures/Figure_3/generate_figure_3_abc.R), panels A-C |
| [`Alpha_Diversity.R`](Alpha_Diversity.R) | [`generate_figure_3_dg.R`](../../figures/Figure_3/generate_figure_3_dg.R), panels D-G, and [`generate_figure_s1_ah.R`](../../supplementary_figures/Figure_S1/generate_figure_s1_ah.R), supplementary panels A-H |
| [`Abundancia_boxplots.R`](Abundancia_boxplots.R) | [`generate_figure_3_hk.R`](../../figures/Figure_3/generate_figure_3_hk.R), panels H-K; [`generate_figure_s2_abc.R`](../../supplementary_figures/Figure_S2/generate_figure_s2_abc.R), supplementary panels A-C; [`generate_figure_s3_aq.R`](../../supplementary_figures/Figure_S3/generate_figure_s3_aq.R), panels A-Q; [`generate_figure_s4_a_an.R`](../../supplementary_figures/Figure_S4/generate_figure_s4_a_an.R), panels A-AN; and [`generate_figure_s5_ao.R`](../../supplementary_figures/Figure_S5/generate_figure_s5_ao.R), panels A-O |
| [`Heatmap_Zscore_MestradoAlberto.R`](Heatmap_Zscore_MestradoAlberto.R) | Figure 4 generators [`A`](../../figures/Figure_4/generate_figure_4_a.R), [`B`](../../figures/Figure_4/generate_figure_4_b.R), [`C`](../../figures/Figure_4/generate_figure_4_c.R), and [`D`](../../figures/Figure_4/generate_figure_4_d.R) |
| [`Correlacoes_Metadados.R`](Correlacoes_Metadados.R) | [`generate_figure_5_d.R`](../../figures/Figure_5/generate_figure_5_d.R), panel D |
| [`Picrust_mestrado_Alberto.R`](Picrust_mestrado_Alberto.R) | Figure 6 generators [`A`](../../figures/Figure_6/generate_figure_6_a.R), [`B`](../../figures/Figure_6/generate_figure_6_b.R), and [`C`](../../figures/Figure_6/generate_figure_6_c.R), with shared corrected multi-group logic in [`figure_6_common.R`](../../figures/Figure_6/figure_6_common.R); the [`Figure S7 PCA generator`](../../supplementary_figures/Figure_S7/generate_figure_s7_ad.R); the [`KEGG-pathway KO z-score source-data generator`](../04_functional_prediction/generate_kegg_pathway_ko_z_scores.R) for Figures 7, 8, and S8-S16; the panel-E renderers for [`Figure 7`](../../figures/Figure_7/generate_figure_7_e.R) and [`Figure 8`](../../figures/Figure_8/generate_figure_8_e.R); the shared [`Figure S8-S16 heatmap renderer`](../../supplementary_figures/figure_s8_s16_common.R) with its [`validator`](../../supplementary_figures/validate_figure_s8_s16.R); and the [`Figure S18 generator`](../../supplementary_figures/Figure_S18/generate_figure_s18.R) |
| [`Script_picrust_mehg.txt`](Script_picrust_mehg.txt) | [`picrust2_workflow.md`](../04_functional_prediction/picrust2_workflow.md), exact PICRUSt2 inputs and analysis-ready functional-prediction exports used by the Figure 6 workflow |

`Script_picrust_mehg.txt` documents successive batch 02 attempts in a command log that was overwritten during analysis. The portable workflow identifies the final successful filtered input, corrects the recorded sample/ASV-output counts from the supplied files, and keeps the historical text unchanged so the reconstruction remains auditable.

The Dams-F0 nonparametric block in `Abundancia_boxplots.R` supplies the
analysis and visual provenance for Figure S2. It read the historical
`AbundanciasTaxons_Todos.csv` and `Metadata_mehg.csv` files, selected
`Geracao_Sexo == "F0_Dams"`, converted taxon labels with `make.names()`, used a
four-group Kruskal-Wallis test, and annotated FDR-adjusted Dunn comparisons.
Its three source panels were exported as 600 x 500 TIFFs under the archived
`Alberto/images/Abundancias/NONparametric/Dams_F0_NONparametric/` directory and
later assembled in Photoshop. The portable
[`generate_figure_s2_abc.R`](../../supplementary_figures/Figure_S2/generate_figure_s2_abc.R)
reads the already normalized batch-01 abundance table and canonical manifest,
reproduces those nonparametric statistics and dotted display labels, and emits
only the three individual PNG panels through the shared
[`figure_s2_s5_abundance_common.R`](../../supplementary_figures/figure_s2_s5_abundance_common.R)
workflow. It neither sources the historical script nor recreates the Photoshop
composite. The paired
[`validator`](../../supplementary_figures/Figure_S2/validate_figure_s2_abc.R)
locks the cohort, taxa, exact test results, significant comparisons, and output
dimensions for traceability.

The PND21 and one-year nonparametric blocks in `Abundancia_boxplots.R` also
supply the analysis and visual provenance for Figures S3-S5. Figure S3 selects
female and male F1-PND21 samples from batch 01; Figures S4 and S5 select female
and male F1-PND365 and F2-PND365 samples, respectively, from batch 02. The
historical numerical inputs are byte-identical to the canonical batch-specific
relative-abundance tables, and the portable workflows join them to the
canonical manifests rather than duplicating source data. All four supplementary
generators share the validated
[`figure_s2_s5_abundance_common.R`](../../supplementary_figures/figure_s2_s5_abundance_common.R)
renderer and preserve the legacy four-group Kruskal-Wallis plus six-pair Dunn
test with FDR adjustment. They recreate all 72 individual source panels as PNG
files without sourcing the legacy script or assembling the Photoshop
composites. Every Figure S3 panel uses the repository-wide treatment order
`CN`, `VitA`, `MeHg`, `VitA+MeHg` and the shared `#98F5FF`, `#B3EE3A`,
`#CD2626`, and `#4D4D4D` fills.
Figure S3-S5 female panels are written to `output/females/` and
male panels to `output/male/`, while each final composite remains at its output
root. Their paired validators lock source checksums, cohorts, taxa, exact
nonparametric results, significant comparison signatures, filenames,
sex-specific paths, and dimensions.

`Picrust_mestrado_Alberto.R` is preserved unchanged. The portable Figure 6 workflow uses only its KO-to-KEGG, four-group ALDEx2, six-comparison Bonferroni, pathway-selection, and visual-design provenance. It intentionally omits the exploratory PCA, EC, MetaCyc, GSEA, z-score heatmap, and mercuric-KO sections from the Figure 6 runtime. Its edited `pathway_errorbar3()` was not sourced because the legacy implementation could combine differently filtered/ordered DAA tables and was not safe for four-group fold-change bookkeeping; the corrected implementation aligns all rows by explicit keys and disables the undefined two-group fold-change panel.

The PCA block around lines 658-742 supplies the analytical and visual
provenance for
[`generate_figure_s7_ad.R`](../../supplementary_figures/Figure_S7/generate_figure_s7_ad.R).
The portable script retains `prcomp()` on transposed, untransformed KO scores,
centering without variance scaling, 68% treatment ellipses, the ggplot2
four-group colors, and sex point shapes. It replaces the unavailable
`ggbiplot` dependency with equivalent `ggplot2` geometry, aligns every matrix
column through the canonical manifests, makes the `10F`, `22F`, and `t3`
exclusions explicit, and fixes PCA axis signs to the supplied reference. The
legacy file remains unchanged and is not sourced by the portable workflow.

The KO-level z-score calculation previously embedded near the end of the legacy script is isolated in [`generate_kegg_pathway_ko_z_scores.R`](../04_functional_prediction/generate_kegg_pathway_ko_z_scores.R), which uses the canonical manifest, descriptive paths, explicit validation, and clean CSV output. The legacy `gplots::heatmap.2()` calls around lines 1003-1065 provide the rendering provenance for [`generate_figure_7_e.R`](../../figures/Figure_7/generate_figure_7_e.R) and [`generate_figure_8_e.R`](../../figures/Figure_8/generate_figure_8_e.R). The portable panel-E scripts read the clean F1-PND21/map00310 and map00480 sources, preserve complete-linkage row clustering and the 15-step blue-white-red palette, omit column clustering and the key as in the source panels used during migration, and record each reference panel's legal dendrogram rotation to prevent platform-dependent branch flips among tied standardized row means. The standalone Figure 8 heatmap TIFF used during migration is not retained. The scripts do not source or execute the historical script.

The mercuric-KO block around lines 1079-1143 supplies the KO-selection and
visual-design provenance for
[`generate_figure_s18.R`](../../supplementary_figures/Figure_S18/generate_figure_s18.R).
It queried seven annotations containing "mercuric" and plotted three KOs
(`K00520`, `K08365`, and `K19058`), used
the treatment palette `#98F5FF`, `#B3EE3A`, `#CD2626`, and `#4D4D4D`, and
faceted the boxplots by generation and sex. The historical block is not
independently executable because it reads from an object named `KO` that is
never created in the script. The official generator therefore reads both
canonical PICRUSt2 matrices and manifests directly, uses the four KOs common to
the selected batch 01 and batch 02 cohorts (`K00520`, `K08363`, `K08364`, and
`K08365`), and compares the F0-Dams, F1-PND21, F1-PND365, and F2-PND365
datasets in a 4 x 4 layout. It retains the historical exposure colors, displays
sex as point shape, and adds independently validated one-way treatment ANOVA
and Tukey HSD calculations without sourcing or mutating the legacy file. The
manually retrieved KEGG result is retained under
[`data/reference/kegg/`](../../data/reference/kegg/) with its original contents
and a repository-compliant filename.

Unreferenced exploratory scripts formerly stored directly under `scripts/` were removed during the repository-wide review. Every file retained in this directory is explicitly mapped to a maintained workflow or documented downstream provenance above.

The three preprocessing-related historical scripts were moved here without changing their contents. Their hard-coded user paths, package-installation calls, interactive `View()` calls, TIFF exports, whole-session saves, inconsistent taxonomy inputs, and exploratory downstream blocks are retained only for auditability. The portable replacement uses the DOI-pinned SILVA 138.1 DADA2 FASTAs with the historical `assignTaxonomy()` and `addSpecies()` method, writes repository-standard data products, and keeps raw reads and large references outside version control. See the [`16S preprocessing guide`](../01_sequence_processing/README.md) for the detailed mapping and current validation boundary.

The historical provenance spans R, DADA2, phyloseq, tidyverse/ggplot2, PICRUSt2, ggpicrust2, ALDEx2, SILVA, and KEGG [@r-core; @dada2; @phyloseq; @tidyverse; @ggplot2; @picrust2; @ggpicrust2; @aldex2; @silva; @kegg].

## References

- <a id="ref-dada2"></a>Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods, 13*(7), 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- <a id="ref-picrust2"></a>Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., Huttenhower, C., & Langille, M. G. I. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology, 38*(6), 685–688. [https://doi.org/10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
- <a id="ref-aldex2"></a>Fernandes, A. D., Reid, J. N. S., Macklaim, J. M., McMurrough, T. A., Edgell, D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput sequencing datasets: Characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome, 2*, 15. [https://doi.org/10.1186/2049-2618-2-15](https://doi.org/10.1186/2049-2618-2-15)
- <a id="ref-kegg"></a>Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27–30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- <a id="ref-phyloseq"></a>McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE, 8*(4), e61217. [https://doi.org/10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
- <a id="ref-silva"></a>Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. *Nucleic Acids Research, 41*(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
- <a id="ref-ggpicrust2"></a>Yang, C., Mai, J., Cao, X., Burberry, A., Cominelli, F., & Zhang, L. (2023). ggpicrust2: An R package for PICRUSt2 predicted functional profile analysis and visualization. *Bioinformatics, 39*(8), btad470. [https://doi.org/10.1093/bioinformatics/btad470](https://doi.org/10.1093/bioinformatics/btad470)
- <a id="ref-ggplot2"></a>Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis* (2nd ed.). Springer. [https://doi.org/10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4)
- <a id="ref-tidyverse"></a>Wickham, H., Averick, M., Bryan, J., et al. (2019). Welcome to the Tidyverse. *Journal of Open Source Software, 4*(43), 1686. [https://doi.org/10.21105/joss.01686](https://doi.org/10.21105/joss.01686)
