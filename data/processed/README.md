# Processed microbiome data

Processed sequence products are separated by sequencing batch. Filenames describe their contents while preserving the file contents from the source folders.

- `asv_counts.tsv`: ASV-by-sample count table.
- `asv_taxonomy.tsv`: ASV taxonomy assignments.
- `asv_counts_with_taxonomy.tsv`: combined taxonomy and count table, when available.
- `asv_relative_abundance_with_taxonomy.tsv`: combined taxonomy and relative-abundance table, when available.
- `representative_sequences.fasta`: representative ASV sequences.
- `feature_table.biom`: BIOM-format feature table.
- `read_tracking.csv`: read counts across filtering/denoising stages, when available.

The data are archived without removing excluded or unexpected sample columns. Use the canonical sample manifest in `data/metadata/` to select analysis samples.

## Raw-sequence availability

The underlying raw reads are not stored in this repository. They are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample ranges are `SAMN61406324`–`SAMN61406377` for Dams_F0/F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365/F2-PND365. The records are on hold and will be released publicly upon manuscript publication.
