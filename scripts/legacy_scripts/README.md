# Legacy scripts used for analysis and figure provenance

These historical scripts were used as references when portable analysis and figure workflows were created. Their contents are preserved unchanged for provenance; they may contain historical paths, interactive state assumptions, package-installation calls, overwritten exploratory commands, and obsolete code, so validated workflows do not source or execute them at runtime.

| Legacy reference script | Portable workflow derived from it |
|---|---|
| [`Abundancias_Barras_Mestrado.R`](Abundancias_Barras_Mestrado.R) | [`generate_figure_3_abc.R`](../../figures/Figure_3/generate_figure_3_abc.R), panels A-C |
| [`Alpha_Diversity.R`](Alpha_Diversity.R) | [`generate_figure_3_dg.R`](../../figures/Figure_3/generate_figure_3_dg.R), panels D-G |
| [`Abundancia_boxplots.R`](Abundancia_boxplots.R) | [`generate_figure_3_hk.R`](../../figures/Figure_3/generate_figure_3_hk.R), panels H-K |
| [`Heatmap_Zscore_MestradoAlberto.R`](Heatmap_Zscore_MestradoAlberto.R) | Figure 4 generators [`A`](../../figures/Figure_4/generate_figure_4_a.R), [`B`](../../figures/Figure_4/generate_figure_4_b.R), [`C`](../../figures/Figure_4/generate_figure_4_c.R), and [`D`](../../figures/Figure_4/generate_figure_4_d.R) |
| [`Correlacoes_Metadados.R`](Correlacoes_Metadados.R) | [`generate_figure_5_d.R`](../../figures/Figure_5/generate_figure_5_d.R), panel D |
| [`Picrust_mestrado_Alberto.R`](Picrust_mestrado_Alberto.R) | Figure 6 generators [`A`](../../figures/Figure_6/generate_figure_6_a.R), [`B`](../../figures/Figure_6/generate_figure_6_b.R), and [`C`](../../figures/Figure_6/generate_figure_6_c.R), with shared corrected multi-group logic in [`figure_6_common.R`](../../figures/Figure_6/figure_6_common.R) |
| [`Script_picrust_mehg.txt`](Script_picrust_mehg.txt) | [`picrust2_workflow.md`](../04_functional_prediction/picrust2_workflow.md), exact PICRUSt2 inputs and analysis-ready functional-prediction exports used by the Figure 6 workflow |

`Script_picrust_mehg.txt` documents successive batch 02 attempts in a command log that was overwritten during analysis. The portable workflow identifies the final successful filtered input, corrects the recorded sample/ASV-output counts from the supplied files, and keeps the historical text unchanged so the reconstruction remains auditable.

`Picrust_mestrado_Alberto.R` is preserved unchanged. The portable Figure 6 workflow uses only its KO-to-KEGG, four-group ALDEx2, six-comparison Bonferroni, pathway-selection, and visual-design provenance. It intentionally omits the exploratory PCA, EC, MetaCyc, GSEA, z-score heatmap, and mercuric-KO sections. Its edited `pathway_errorbar3()` was not sourced because the legacy implementation could combine differently filtered/ordered DAA tables and was not safe for four-group fold-change bookkeeping; the corrected implementation aligns all rows by explicit keys and disables the undefined two-group fold-change panel.

Legacy scripts that have not yet been used for a portable workflow remain directly under [`scripts/`](../). This separation makes it possible to identify migrated provenance files without changing or prematurely curating the unused scripts.
