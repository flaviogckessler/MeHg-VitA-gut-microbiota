# Legacy KO z-score validation checkpoints

These 19 CSV files are the unchanged exports copied from the historical
`calculate_ZScore_KO2kegg()` workflow in
[`Picrust_mestrado_Alberto.R`](../../../scripts/legacy_scripts/Picrust_mestrado_Alberto.R).
They preserve the KO-level treatment z-scores used for the Figure 7, Figure 8,
and Figure S8-S16 heatmaps.

They are validation checkpoints, not primary figure inputs. The portable
[`generate_kegg_pathway_ko_z_scores.R`](../../../scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R)
workflow regenerates cleanly named source tables under each batch directory in
[`results/functional_prediction/`](../../functional_prediction/) using the
canonical PICRUSt2 KO matrices and without reading these checkpoints. The clean
tables omit the legacy row-number column and rename `ID` to `ko_id`; KO order,
missing-value pattern, and z-score values are otherwise preserved.

The matching
[`validate_kegg_pathway_ko_z_scores.R`](../../../scripts/04_functional_prediction/validate_kegg_pathway_ko_z_scores.R)
script compares all 19 regenerated tables with these checkpoints at a numeric
tolerance of `1e-12`. The initial migration comparison covered 229 pathway-KO
rows and had a maximum absolute difference of `9.992007e-15`, attributable only
to floating-point summation order.

The legacy filename suffix `1year` corresponds to PND365. The double underscore
in filenames such as `ko00190__F1_PND21.csv` is retained here because these files
are unchanged provenance artifacts.

The Figure S16 checkpoint is `ko00311__F2_1year.csv`, corresponding to the
F2-PND365 cohort used for that figure.

## References

These checkpoints derive from R-based ggpicrust2 and KEGG processing [@r-core; @ggpicrust2; @kegg]. Complete entries are in [`bibliography.bib`](../../../bibliography.bib).
