# Data organization

The microbiome data are organized by sequencing batch and processing level.

## Sequencing batches

| Batch directory | Biological material | Time point | Sequencing platform |
|---|---|---|---|
| `batch01_iontorrent_dams-f0_f1-pnd21` | F0 dams and F1 offspring | F1 at PND21; dam age not encoded | Ion Torrent |
| `batch02_illumina-miseq_f1-f2_pnd365` | F1 and F2 offspring | PND365 | Illumina MiSeq |

## Directories

- `metadata/` contains preserved source metadata, standardized per-batch sample manifests, a combined sample manifest, SRA deposition information, and validation results.
- `processed/` contains ASV count tables, taxonomy, representative sequences, BIOM files, relative abundance data, and read tracking.
- `reference/` contains small, versioned reference objects required to run portable analyses, including the repository-local ggpicrust2 KEGG-to-KO mapping used by the KO z-score generator and the preserved KEGG REST query result used to select the mercuric-resistance KOs in Figure S18.
- `external/` documents the ignored local location for large external reference databases that must not be committed. The required SILVA 138.1 taxonomy and species FASTAs are downloaded from the versioned [Zenodo record](https://doi.org/10.5281/zenodo.4587955); filenames, checksums, licensing, and download instructions are recorded in [`external/README.md`](external/README.md).

Raw FASTQ files are deliberately not stored in `New_repo`, and there is no local `data/raw/` directory.

Analysis scripts should join data to `metadata/sample_manifest_all_batches.csv` by `sample_id` and should honor `include_primary_analysis`. Samples `10F` and `22F` are retained in the archived data but flagged for exclusion because of low sequencing quality.

The batch 02 processed count table contains `t3`, although `t3` is absent from the provided sample list and metadata. It is preserved in the source-derived table and reported as unexpected in `metadata/data_validation_summary.csv`; analyses should not include it unless valid metadata are supplied.

## Raw-sequence availability

Raw sequencing data generated during this study have been deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample accessions are `SAMN61406324`–`SAMN61406377` for Dams_F0 and F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365 and F2-PND365. The records are currently on hold and will be released publicly upon publication of the manuscript.

The accession ranges are also recorded in `metadata/sra_data_availability.csv`. No one-to-one sample-to-BioSample mapping is asserted because that mapping was not supplied.

For sequence reprocessing after public release, download the FASTQ files outside the repository and follow the [`scripts/01_sequence_processing/` guide](../scripts/01_sequence_processing/README.md). Set `DADA2_BATCH_ID` and `RAW_FASTQ_DIR` before running `scripts/01_sequence_processing/preprocess_16s.R`. The script does not assume or create a repository-local raw-data directory, and it writes reviewed reruns below the selected processed batch without overwriting archived canonical products.

The historical [`scripts/00_setup/organize_data.ps1`](../scripts/00_setup/organize_data.ps1) utility documents the initial curation process. Portable analyses should use the canonical files already organized under `data/` and `results/` rather than rebuilding them through the historical import workflow.

Relevant records cover DADA2, SILVA 138.1, the versioned DADA2 training dataset, PICRUSt2, KEGG, and the Sequence Read Archive [@dada2; @silva; @silva-dada2-138-1; @picrust2; @kegg; @sra].

## References

- <a id="ref-dada2"></a>Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods, 13*(7), 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- <a id="ref-picrust2"></a>Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., Huttenhower, C., & Langille, M. G. I. (2020). PICRUSt2 for prediction of metagenome functions. *Nature Biotechnology, 38*(6), 685–688. [https://doi.org/10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
- <a id="ref-kegg"></a>Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research, 28*(1), 27–30. [https://doi.org/10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27)
- <a id="ref-sra"></a>Leinonen, R., Sugawara, H., Shumway, M., & International Nucleotide Sequence Database Collaboration. (2011). The Sequence Read Archive. *Nucleic Acids Research, 39*(Database issue), D19–D21. [https://doi.org/10.1093/nar/gkq1019](https://doi.org/10.1093/nar/gkq1019)
- <a id="ref-silva-dada2-138-1"></a>McLaren, M. R., & Callahan, B. J. (2021). *Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2* (Version 1) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.4587955](https://doi.org/10.5281/zenodo.4587955)
- <a id="ref-silva"></a>Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. *Nucleic Acids Research, 41*(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
