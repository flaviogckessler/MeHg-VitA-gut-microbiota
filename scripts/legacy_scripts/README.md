# Legacy scripts used for figure provenance

These historical scripts were used as references when the portable Figure 3, Figure 4, and Figure 5D workflows were created. Their contents are preserved unchanged for provenance; they may contain historical paths, interactive state assumptions, package-installation calls, and obsolete exploratory code, so the validated figure generators do not source or execute them at runtime.

| Legacy reference script | Portable workflow derived from it |
|---|---|
| [`Abundancias_Barras_Mestrado.R`](Abundancias_Barras_Mestrado.R) | [`generate_figure_3_abc.R`](../../figures/Figure_3/generate_figure_3_abc.R), panels A-C |
| [`Alpha_Diversity.R`](Alpha_Diversity.R) | [`generate_figure_3_dg.R`](../../figures/Figure_3/generate_figure_3_dg.R), panels D-G |
| [`Abundancia_boxplots.R`](Abundancia_boxplots.R) | [`generate_figure_3_hk.R`](../../figures/Figure_3/generate_figure_3_hk.R), panels H-K |
| [`Heatmap_Zscore_MestradoAlberto.R`](Heatmap_Zscore_MestradoAlberto.R) | Figure 4 generators [`A`](../../figures/Figure_4/generate_figure_4_a.R), [`B`](../../figures/Figure_4/generate_figure_4_b.R), [`C`](../../figures/Figure_4/generate_figure_4_c.R), and [`D`](../../figures/Figure_4/generate_figure_4_d.R) |
| [`Correlacoes_Metadados.R`](Correlacoes_Metadados.R) | [`generate_figure_5_d.R`](../../figures/Figure_5/generate_figure_5_d.R), panel D |

Legacy scripts that have not yet been used for a portable workflow remain directly under [`scripts/`](../). This separation makes it possible to identify migrated provenance files without changing or prematurely curating the unused scripts.
