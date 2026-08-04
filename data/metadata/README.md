# Metadata files

`sample_manifest_all_batches.csv` is the canonical sample-level metadata table for analysis. It standardizes batch, platform, generation, life stage, collection day, treatment, sex, analysis inclusion, SRA deposition provenance, release status, and source metadata across both sequencing batches.

Each batch directory contains:

- an unchanged copy of the provided source metadata;
- a standardized `sample_manifest.csv`;
- any additional legacy metadata available for that batch.

`data_validation_summary.csv` compares the expected sample lists, metadata rows, and processed ASV-table columns. Empty missing/unexpected fields mean that the identifiers agree.

The `legacy_observation` column preserves the unexplained `obs` field from the batch 02 source metadata without assigning it a new interpretation.

For the 12 combined-treatment samples in batch 01, the designated source metadata uses `VitA+MeHg` while the legacy tab-separated metadata uses `VitA|MeHg`. The canonical manifest follows the designated `Metadata_mehg.csv` spelling; the unchanged legacy file is retained for provenance.

## Raw-sequence availability

Raw reads are not stored in this repository. They have been deposited in the NCBI Sequence Read Archive under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663):

- Dams_F0 and F1-PND21: `SAMN61406324`–`SAMN61406377`.
- F1-PND365 and F2-PND365: `SAMN61407672`–`SAMN61407719`.

The records are currently on hold and will be made public upon publication of the manuscript. `sra_data_availability.csv` provides these batch-level ranges in machine-readable form. It does not infer individual sample-to-BioSample assignments.
