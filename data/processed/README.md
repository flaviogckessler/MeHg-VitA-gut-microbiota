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

## FASTQ reprocessing outputs

The canonical [`preprocess_16s.R`](../../scripts/01_sequence_processing/preprocess_16s.R) workflow writes new results to `<batch>/dada2_rerun/` by default. This review location uses the same descriptive filenames listed above but protects the supplied archived products from accidental replacement. It also records `preprocessing_parameters.tsv` and `preprocessing_session_info.txt`. Promote rerun files to the batch root only after FASTQ availability, sample-identity validation, read-retention review, and regression comparison with the archived products.

The workflow reproduces the batch-specific archived schemas rather than creating a synthetic common schema. Batch 01 produces the combined count/taxonomy and renamed-taxonomy files. Batch 02 produces the filtered relative-abundance/taxonomy file. New products are compared automatically with this directory by a local-only preprocessing validator retained outside version control; a rerun fails if the canonical contents differ. BIOM comparison is semantic because generated JSON timestamps are not scientifically meaningful.

Exact batch 02 reproduction retains the historical `t3` column during preprocessing, even though `t3` is absent from the supplied metadata. This is an output-compatibility decision, not permission to analyze the sample. Downstream workflows must exclude `t3` unless valid metadata are supplied.

The batch 02 relative-abundance table is a filtered derivative of the 11,704-ASV count table. It contains 2,739 bacterial ASVs classified through Family and passing the documented prevalence/abundance thresholds. `SV` numbers refer to positions in the unfiltered table, and relative abundances are calculated after filtering. See the [preprocessing guide](../../scripts/01_sequence_processing/README.md) for the exact rule and the preserved historical genus-label convention.

## PICRUSt2-specific inputs

The `picrust2_input_*` files are preserved separately because they record the exact pipeline inputs and must not be reconstructed by silently substituting the broader canonical files:

- Batch 01 used 1,684 ASVs. Its input is not an exact subset of the 1,701-ASV canonical table: 22 sequence identifiers are absent from the canonical file and 1,065 shared rows have different count vectors.
- Batch 02 used a filtered set of 2,739 of the 11,704 canonical sequences. Its `SV00001`-style FASTA identifiers link one-to-one to the count-table identifiers; resolving those identifiers through the FASTA confirms that all 2,739 sequences and count vectors match the canonical table.

Both PICRUSt2 count tables retain all columns used in the historical runs. This includes excluded batch 01 samples `10F` and `22F` and the undocumented batch 02 column `t3`; downstream analyses must apply the canonical manifest. The intermediate HDF5 BIOM tables are reproducibly generated from these TSV files and are not duplicated in the repository. See [`scripts/04_functional_prediction/picrust2_workflow.md`](../../scripts/04_functional_prediction/picrust2_workflow.md).

## Raw-sequence availability

The underlying raw reads are not stored in this repository. They are deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample ranges are `SAMN61406324`–`SAMN61406377` for Dams_F0/F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365/F2-PND365. The records are on hold and will be released publicly upon manuscript publication.

The processed products and interchange formats are associated with DADA2, phyloseq, the `biomformat` interface, BIOM, PICRUSt2, and the Sequence Read Archive [@dada2; @phyloseq; @biomformat; @biom; @picrust2; @sra].

## References

- <a id="ref-dada2"></a>Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods, 13*(7), 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- <a id="ref-picrust2"></a>Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., Huttenhower, C., & Langille, M. G. I. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology, 38*(6), 685–688. [https://doi.org/10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
- <a id="ref-biom"></a>McDonald, D., Clemente, J. C., Kuczynski, J., Rideout, J. R., Stombaugh, J., Wendel, D., Wilke, A., Huse, S., Hufnagle, J., Meyer, F., Knight, R., & Caporaso, J. G. (2012). The Biological Observation Matrix (BIOM) format or: How I learned to stop worrying and love the ome-ome. *GigaScience, 1*(1), 7. [https://doi.org/10.1186/2047-217X-1-7](https://doi.org/10.1186/2047-217X-1-7)
- <a id="ref-biomformat"></a>McMurdie, P. J., & Paulson, J. N. (2026). *biomformat: An interface package for the BIOM file format* (Version 1.40.0) [R/Bioconductor package]. [https://doi.org/10.18129/B9.bioc.biomformat](https://doi.org/10.18129/B9.bioc.biomformat)
- <a id="ref-phyloseq"></a>McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE, 8*(4), e61217. [https://doi.org/10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
- <a id="ref-sra"></a>Leinonen, R., Sugawara, H., Shumway, M., & International Nucleotide Sequence Database Collaboration. (2011). The Sequence Read Archive. *Nucleic Acids Research, 39*(Database issue), D19–D21. [https://doi.org/10.1093/nar/gkq1019](https://doi.org/10.1093/nar/gkq1019)
