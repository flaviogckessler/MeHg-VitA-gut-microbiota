# Results

- `alpha_diversity/` contains batch-specific alpha-diversity estimates and their source-associated metadata. Figure workflows join these estimates to the canonical manifests in `data/metadata/`.
- `taxon_relative_abundance/` contains batch-specific processed taxon relative-abundance matrices used by Figure 3 panels H-K.
- `functional_prediction/` contains legacy PICRUSt exports, separated by sequencing batch and renamed by content.
- `legacy_summary_tables/` contains the cross-batch Z-score table found in the original `Data_metadata` folder.

The `kegg_ortholog_abundance_legacy_ko2.tsv` filenames deliberately retain the `KO2` distinction because the precise transformation represented by `KO2` is not documented in the source files.

## Raw-sequence availability

The raw reads underlying these results are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample accessions span `SAMN61406324`–`SAMN61406377` for Dams_F0/F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365/F2-PND365. They are currently on hold and will become public upon publication of the manuscript.
