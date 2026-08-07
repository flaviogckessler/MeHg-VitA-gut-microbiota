# Processed microbiome data

Processed sequence products are separated by sequencing batch. Filenames describe their contents while preserving the file contents from the source folders.

- `asv_counts.tsv`: ASV-by-sample count table.
- `asv_taxonomy.tsv`: ASV taxonomy assignments.
- `asv_counts_with_taxonomy.tsv`: combined taxonomy and count table, when available.
- `asv_relative_abundance_with_taxonomy.tsv`: combined taxonomy and relative-abundance table, when available.
- `representative_sequences.fasta`: representative ASV sequences.
- `feature_table.biom`: BIOM-format feature table.
- `read_tracking.csv`: read counts across filtering/denoising stages, when available.
- `picrust2_input_asv_counts.tsv`: exact ASV count table supplied to the retained PICRUSt2 run.
- `picrust2_input_representative_sequences.fasta`: exact representative sequences supplied to the retained PICRUSt2 run.

The data are archived without removing excluded or unexpected sample columns. Use the canonical sample manifest in `data/metadata/` to select analysis samples.

## PICRUSt2-specific inputs

The `picrust2_input_*` files are preserved separately because they record the exact pipeline inputs and must not be reconstructed by silently substituting the broader canonical files:

- Batch 01 used 1,684 ASVs. Its input is not an exact subset of the 1,701-ASV canonical table: 22 sequence identifiers are absent from the canonical file and 1,065 shared rows have different count vectors.
- Batch 02 used a filtered set of 2,739 of the 11,704 canonical sequences. Its `SV00001`-style FASTA identifiers link one-to-one to the count-table identifiers; resolving those identifiers through the FASTA confirms that all 2,739 sequences and count vectors match the canonical table.

Both PICRUSt2 count tables retain all columns used in the historical runs. This includes excluded batch 01 samples `10F` and `22F` and the undocumented batch 02 column `t3`; downstream analyses must apply the canonical manifest. The intermediate HDF5 BIOM tables are reproducibly generated from these TSV files and are not duplicated in the repository. See [`scripts/04_functional_prediction/picrust2_workflow.md`](../../scripts/04_functional_prediction/picrust2_workflow.md).

## Raw-sequence availability

The underlying raw reads are not stored in this repository. They are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample ranges are `SAMN61406324`–`SAMN61406377` for Dams_F0/F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365/F2-PND365. The records are on hold and will be released publicly upon manuscript publication.
