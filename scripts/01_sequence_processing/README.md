# 16S sequence preprocessing

[`preprocess_16s.R`](preprocess_16s.R) is the canonical, linear FASTQ-to-ASV workflow for both sequencing batches. It replaces several non-portable historical scripts while preserving their batch-specific DADA2 settings.

## Raw/source sequencing data

The source sequencing data (raw FASTQ reads) have been deposited in the NCBI Sequence Read Archive (SRA) under [BioProject PRJNA1489663](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1489663). They are not stored in this repository. The records are currently on hold and are intended for public release with the manuscript.

The supplied deposition metadata do not include a one-to-one mapping from repository `sample_id` values to BioSample or SRA run accessions. After public release, obtain that mapping before renaming or processing downloaded files. The script requires FASTQ-derived identifiers to match the canonical batch manifest so an accession cannot silently be assigned to the wrong sample.

## Provenance and decisions

| Batch | Historical source | Preserved settings |
|---|---|---|
| `batch01_iontorrent_dams-f0_f1-pnd21` | `scripts/legacy_scripts/dada2_phyloseq_batch01_iontorrent.R` | Single-end; `truncLen=275`, `trimLeft=20`, `maxEE=2`, pseudo-pooling, and the DADA2 Ion Torrent settings `HOMOPOLYMER_GAP_PENALTY=-1`, `BAND_SIZE=32`. |
| `batch02_illumina-miseq_f1-f2_pnd365` | `scripts/legacy_scripts/dada2_batch02_illumina.R` | Paired-end; `truncLen=c(149,70)`, `trimLeft=c(20,0)`, `maxEE=c(2,2)`, pseudo-pooling, paired-read merging, and consensus chimera removal. |

The historical `phyloseq_exploratory.R` starts from an already constructed `phyloseq` object and performs downstream filtering, abundance summaries, diversity analyses, plotting, and a random-tree experiment. It is not sourced by the preprocessing workflow.

Taxonomy uses the versioned [SILVA 138.1 prokaryotic DADA2 reference](https://doi.org/10.5281/zenodo.4587955). This matches the historical `assignTaxonomy()` followed by `addSpecies()` method. The former supplied `SILVA_SSU_r138_2019.RData` was a different DECIPHER-format classifier and is no longer a workflow input.

## Reference download

Large reference files remain outside version control. From the repository root, download and verify the two required Zenodo files with R:

```r
reference_dir <- file.path("data", "external", "silva")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

download.file(
  "https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1",
  file.path(reference_dir, "silva_nr99_v138.1_train_set.fa.gz"),
  mode = "wb"
)
download.file(
  "https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1",
  file.path(reference_dir, "silva_species_assignment_v138.1.fa.gz"),
  mode = "wb"
)

observed <- unname(tools::md5sum(c(
  file.path(reference_dir, "silva_nr99_v138.1_train_set.fa.gz"),
  file.path(reference_dir, "silva_species_assignment_v138.1.fa.gz")
)))
stopifnot(identical(observed, c(
  "6b41db7139834c71171f8ce5b5918fc6",
  "f21c2d97c79ff07c17949a9622371a4c"
)))
```

These MD5 values are published with Zenodo record `4587955`. See [`data/external/README.md`](../../data/external/README.md) for file sizes, licensing, alternative locations, and a documented caveat affecting a small number of taxonomy labels.

## Inputs and filenames

Install `dada2` and `biomformat` from Bioconductor before a full run. Keep raw reads outside the repository. Set `SILVA_TAXONOMY_FASTA` and `SILVA_SPECIES_FASTA` only when the downloaded references are stored somewhere other than `data/external/silva/`.

FASTQ files must be demultiplexed and named with canonical sample identifiers:

- batch 01: `<sample_id>.IonXpress<optional text>.fastq.gz`; the `.IonXpress...` portion is optional, and `.fastq`, `.fq`, and gzip variants are accepted;
- batch 02: `<sample_id>_R1_001.fastq.gz` and `<sample_id>_R2_001.fastq.gz`; `.fastq`, `.fq`, and gzip variants are accepted.

Every sample in the selected batch manifest must be present exactly once. The undocumented historical batch 02 column `t3` is not accepted unless valid canonical metadata are first added.

## Configuration validation without FASTQ data

Validation-only mode parses the script and checks the selected canonical manifest. If the SILVA FASTAs have been downloaded, it also verifies both files against the Zenodo MD5 values. FASTQ files and Bioconductor packages are not required.

```powershell
$env:DADA2_BATCH_ID = "batch01_iontorrent_dams-f0_f1-pnd21"
$env:DADA2_VALIDATE_ONLY = "true"
Rscript scripts/01_sequence_processing/preprocess_16s.R
Remove-Item Env:DADA2_VALIDATE_ONLY
```

Repeat with `batch02_illumina-miseq_f1-f2_pnd365` to validate the second manifest.

## Full run after SRA release

```powershell
$env:DADA2_BATCH_ID = "batch01_iontorrent_dams-f0_f1-pnd21"
$env:RAW_FASTQ_DIR = "D:\path\to\demultiplexed_fastq"
Rscript scripts/01_sequence_processing/preprocess_16s.R
```

Change `DADA2_BATCH_ID` for batch 02. Outputs go to `data/processed/<batch_id>/dada2_rerun/` by default, protecting the archived canonical products from accidental replacement. Set `DADA2_OUTPUT_DIR` for another review location.

The script writes `asv_counts.tsv`, `asv_taxonomy.tsv`, `asv_counts_with_taxonomy.tsv`, `asv_relative_abundance_with_taxonomy.tsv`, `representative_sequences.fasta`, `feature_table.biom`, `read_tracking.csv`, `preprocessing_parameters.tsv`, and `preprocessing_session_info.txt`. The parameters file records the Zenodo DOI and both reference checksums.

## Current validation boundary

The script has been syntax-checked in R and validation-only mode has been run against both manifests and the Zenodo reference configuration. The archived products confirm the expected final sequence lengths: all 1,701 batch 01 ASVs are 255 nt, and all 11,704 batch 02 ASVs are 150 nt.

The denoising and taxonomy run cannot be tested until the SRA FASTQ files and accession-to-sample mapping are public. At that time, review quality profiles, read retention, paired-read overlap, ASV counts, taxonomy, and checksums before replacing archived products. Change trimming values only if the released reads demonstrate that the deposited data differ from the historical inputs.
