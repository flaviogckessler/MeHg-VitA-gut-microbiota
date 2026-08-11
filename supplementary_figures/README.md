# Supplementary figures

## Figures S8-S16

- Their KO-level heatmap source tables are generated under the batch-specific [`results/functional_prediction/`](../results/functional_prediction/) directories by [`generate_kegg_pathway_ko_z_scores.R`](../scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R).
- [`validate_kegg_pathway_ko_z_scores.R`](../scripts/04_functional_prediction/validate_kegg_pathway_ko_z_scores.R) regenerates and checks the complete 19-table source set used across Figures 7, 8, and S8-S16 against the preserved legacy CSV exports.
- As in main Figures 7 and 8, the KEGG pathway portions were prepared manually: one `https://www.kegg.jp/pathway/mapxxxxx` export per exposure group (`CN`, `VitA`, `MeHg`, and `VitA+MeHg`), with pathway boxes colored from the corresponding KO z-score column using the common blue-white-red scale. Here, `xxxxx` is the pathway code represented in the individual figure.
- The supplementary heatmap rendering workflows have not yet been migrated.

## Figure S17

- Source: [`Figure_S17_SourceData.pzfx`](Figure_S17_SourceData.pzfx).
- Generation: the full figure was generated in GraphPad Prism and exported as [`Figure_S17.png`](Figure_S17.png).
