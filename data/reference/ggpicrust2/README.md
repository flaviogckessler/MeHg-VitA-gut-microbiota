# ggpicrust2 KEGG-to-KO reference

`kegg_reference.RData` is a repository-local copy of the offline
KEGG-pathway-to-KO membership reference distributed with `ggpicrust2` 2.0.0.
It is loaded directly by
[`generate_kegg_pathway_ko_z_scores.R`](../../../scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R),
so generating the heatmap source tables does not depend on the user's installed
`ggpicrust2` version or library location.

Provenance:

- package: `ggpicrust2` 2.0.0;
- package-relative source: `extdata/kegg_reference.RData`;
- stored object: `ko_to_kegg_reference`;
- expected dimensions: 306 pathways by 326 columns;
- file size: 62,059 bytes;
- MD5: `67728b7ea0ca954b767f77406e7a44b9`.

The file is copied unchanged. The generator verifies its checksum, object name,
dimensions, pathway IDs, and KO identifiers before using it. The historical CSV
exports under `results/legacy_summary_tables/` are not generation inputs.

This offline mapping is attributed to ggpicrust2 and KEGG [@ggpicrust2; @kegg].

## References

- <a id="ref-kegg"></a>Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27–30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- <a id="ref-ggpicrust2"></a>Yang, C., Mai, J., Cao, X., Burberry, A., Cominelli, F., & Zhang, L. (2023). ggpicrust2: An R package for PICRUSt2 predicted functional profile analysis and visualization. *Bioinformatics, 39*(8), btad470. [https://doi.org/10.1093/bioinformatics/btad470](https://doi.org/10.1093/bioinformatics/btad470)
