# Legacy scripts used for analysis and figure provenance

These historical scripts were used as references when portable analysis and figure workflows were created. Their contents are preserved unchanged for provenance; they may contain historical paths, interactive state assumptions, package-installation calls, overwritten exploratory commands, and obsolete code, so validated workflows do not source or execute them at runtime.

| Legacy reference script | Portable workflow derived from it |
|---|---|
| [`Abundancias_Barras_Mestrado.R`](Abundancias_Barras_Mestrado.R) | [`generate_figure_3_abc.R`](../../figures/Figure_3/generate_figure_3_abc.R), panels A-C |
| [`Alpha_Diversity.R`](Alpha_Diversity.R) | [`generate_figure_3_dg.R`](../../figures/Figure_3/generate_figure_3_dg.R), panels D-G, and [`generate_figure_s1_ah.R`](../../supplementary_figures/Figure_S1/generate_figure_s1_ah.R), supplementary panels A-H |
| [`Abundancia_boxplots.R`](Abundancia_boxplots.R) | [`generate_figure_3_hk.R`](../../figures/Figure_3/generate_figure_3_hk.R), panels H-K |
| [`Heatmap_Zscore_MestradoAlberto.R`](Heatmap_Zscore_MestradoAlberto.R) | Figure 4 generators [`A`](../../figures/Figure_4/generate_figure_4_a.R), [`B`](../../figures/Figure_4/generate_figure_4_b.R), [`C`](../../figures/Figure_4/generate_figure_4_c.R), and [`D`](../../figures/Figure_4/generate_figure_4_d.R) |
| [`Correlacoes_Metadados.R`](Correlacoes_Metadados.R) | [`generate_figure_5_d.R`](../../figures/Figure_5/generate_figure_5_d.R), panel D |
| [`Picrust_mestrado_Alberto.R`](Picrust_mestrado_Alberto.R) | Figure 6 generators [`A`](../../figures/Figure_6/generate_figure_6_a.R), [`B`](../../figures/Figure_6/generate_figure_6_b.R), and [`C`](../../figures/Figure_6/generate_figure_6_c.R), with shared corrected multi-group logic in [`figure_6_common.R`](../../figures/Figure_6/figure_6_common.R); the [`KEGG-pathway KO z-score source-data generator`](../04_functional_prediction/generate_kegg_pathway_ko_z_scores.R) for Figures 7, 8, and S8-S16; the panel-E renderers for [`Figure 7`](../../figures/Figure_7/generate_figure_7_e.R) and [`Figure 8`](../../figures/Figure_8/generate_figure_8_e.R); the shared [`Figure S8-S16 heatmap renderer`](../../supplementary_figures/figure_s8_s16_common.R) with its [`validator`](../../supplementary_figures/validate_figure_s8_s16.R); and the [`Figure S18 generator`](../../supplementary_figures/Figure_S18/generate_figure_s18.R) |
| [`Script_picrust_mehg.txt`](Script_picrust_mehg.txt) | [`picrust2_workflow.md`](../04_functional_prediction/picrust2_workflow.md), exact PICRUSt2 inputs and analysis-ready functional-prediction exports used by the Figure 6 workflow |

`Script_picrust_mehg.txt` documents successive batch 02 attempts in a command log that was overwritten during analysis. The portable workflow identifies the final successful filtered input, corrects the recorded sample/ASV-output counts from the supplied files, and keeps the historical text unchanged so the reconstruction remains auditable.

`Picrust_mestrado_Alberto.R` is preserved unchanged. The portable Figure 6 workflow uses only its KO-to-KEGG, four-group ALDEx2, six-comparison Bonferroni, pathway-selection, and visual-design provenance. It intentionally omits the exploratory PCA, EC, MetaCyc, GSEA, z-score heatmap, and mercuric-KO sections from the Figure 6 runtime. Its edited `pathway_errorbar3()` was not sourced because the legacy implementation could combine differently filtered/ordered DAA tables and was not safe for four-group fold-change bookkeeping; the corrected implementation aligns all rows by explicit keys and disables the undefined two-group fold-change panel.

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

Legacy scripts that have not yet been used for a portable workflow remain directly under [`scripts/`](../). This separation makes it possible to identify migrated provenance files without changing or prematurely curating the unused scripts.
