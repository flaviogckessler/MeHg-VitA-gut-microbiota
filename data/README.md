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
- `external/` is reserved for large external reference databases that should not be committed. No large external reference database is currently stored in the repository.

Raw FASTQ files are deliberately not stored in `New_repo`, and there is no local `data/raw/` directory.

Analysis scripts should join data to `metadata/sample_manifest_all_batches.csv` by `sample_id` and should honor `include_primary_analysis`. Samples `10F` and `22F` are retained in the archived data but flagged for exclusion because of low sequencing quality.

The batch 02 processed count table contains `t3`, although `t3` is absent from the provided sample list and metadata. It is preserved in the source-derived table and reported as unexpected in `metadata/data_validation_summary.csv`; analyses should not include it unless valid metadata are supplied.

## Raw-sequence availability

Raw sequencing data generated during this study have been deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). BioSample accessions are `SAMN61406324`–`SAMN61406377` for Dams_F0 and F1-PND21 and `SAMN61407672`–`SAMN61407719` for F1-PND365 and F2-PND365. The records are currently on hold and will be released publicly upon publication of the manuscript.

The accession ranges are also recorded in `metadata/sra_data_availability.csv`. No one-to-one sample-to-BioSample mapping is asserted because that mapping was not supplied.

For sequence reprocessing after public release, download the FASTQ files outside the repository and set the `RAW_FASTQ_DIR` environment variable before running `scripts/dada2.R`. The script no longer assumes or creates a repository-local raw-data directory.

The historical [`scripts/00_setup/organize_data.ps1`](../scripts/00_setup/organize_data.ps1) utility documents the initial curation process. Portable analyses should use the canonical files already organized under `data/` and `results/` rather than rebuilding them through the historical import workflow.
