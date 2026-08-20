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

The historical `phyloseq_exploratory.R` starts from an already constructed `phyloseq` object and performs downstream filtering, abundance summaries, diversity analyses, plotting, and a random-tree experiment. It is not sourced by the preprocessing workflow. Its documented ASV filter is reproduced only for the batch 02 relative-abundance export so that this derived file matches the archived processed data.

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

Every batch 01 sample in the selected manifest must be present exactly once. Exact reproduction of the archived batch 02 products requires the 48 manifest samples plus the historical `t3` sample, because `t3` is present in every archived batch 02 count vector. The script therefore accepts and requires `t3` for preprocessing compatibility. Its metadata remain unknown, so downstream analyses must continue to exclude it unless valid metadata are supplied.

## Configuration validation without FASTQ data

Validation-only mode parses the script, checks the selected canonical manifest, and audits the internal contract of the archived processed files. The audit verifies ASV/sample dimensions, sequence lengths and order, taxonomy linkage, and batch-specific derived-table rules. If the SILVA FASTAs have been downloaded, it also verifies both files against the Zenodo MD5 values. FASTQ files and Bioconductor packages are not required.

```bash
export DADA2_BATCH_ID="batch01_iontorrent_dams-f0_f1-pnd21"
export DADA2_VALIDATE_ONLY="true"
Rscript scripts/01_sequence_processing/preprocess_16s.R
unset DADA2_VALIDATE_ONLY
```

Repeat with `batch02_illumina-miseq_f1-f2_pnd365` to validate the second manifest.

## Full run after SRA release

```bash
export DADA2_BATCH_ID="batch01_iontorrent_dams-f0_f1-pnd21"
export RAW_FASTQ_DIR="/path/to/demultiplexed_fastq"
Rscript scripts/01_sequence_processing/preprocess_16s.R
```

Change `DADA2_BATCH_ID` for batch 02. Outputs go to `data/processed/<batch_id>/dada2_rerun/` by default, protecting the archived canonical products from accidental replacement. Set `DADA2_OUTPUT_DIR` for another review location.

The canonical output contract is batch-specific:

| Batch | Files that must match the archived products |
|---|---|
| Batch 01 | `asv_counts.tsv`, `asv_taxonomy.tsv`, `asv_counts_with_taxonomy.tsv`, `asv_taxonomy_renamed_ids.tsv`, `representative_sequences.fasta`, `feature_table.biom`, and `read_tracking.csv` |
| Batch 02 | `asv_counts.tsv`, `asv_taxonomy.tsv`, `asv_relative_abundance_with_taxonomy.tsv`, `representative_sequences.fasta`, and `feature_table.biom` |

Batch 01 retains 1,701 255-nt ASVs and 54 samples. Batch 02 retains 11,704 150-nt ASVs and 49 samples in its count, taxonomy, FASTA, and BIOM products. Its relative-abundance compatibility filter keeps Bacteria with assigned Phylum, Class, and Family, prevalence greater than 6 samples, total count greater than 2, and maximum sample count greater than 1. Mitochondrial-family assignments are retained intentionally because the archived tables used by downstream analyses include them. The historical `phyloseq_exploratory.R` condition `Family != "Mitochondrias"` was ineffective against the actual singular value `Mitochondria`; the maintained implementation now expresses the resulting retention policy directly. This yields the same 2,739 ASVs as the archived table, including 17 `Mitochondria` ASVs. Their `SV00001`-style identifiers retain the positions of the unfiltered 11,704-ASV table, including gaps after filtering. Relative abundances are recalculated after filtering, within each sample. The historical display label `Ruminococcaceae Incertae Sedis` is retained in this derived table for 28 genus entries; the full taxonomy export remains unchanged.

The rerun also writes auxiliary `read_tracking.csv` data for batch 02, `preprocessing_parameters.tsv`, and `preprocessing_session_info.txt`. These have no complete archived batch 02 counterparts and are not promoted as historical matches. The parameters file records the Zenodo DOI, checksums, preprocessing settings, and batch 02 compatibility rules.

At the end of a full run, the preprocessing workflow uses a local-only validator to compare generated canonical products with `data/processed/<batch_id>/`. TSV/CSV/FASTA identifiers, order, counts, and taxonomy must agree; batch 02 floating-point relative abundances use a `1e-14` tolerance; BIOM files are compared as feature-by-sample matrices rather than by JSON timestamps. The local validator writes `archived_output_comparison.tsv` and stops with an error if any canonical product differs or a stale cross-batch product is present. It refuses to use the archived batch directory itself as `DADA2_OUTPUT_DIR`, preventing accidental replacement. The standalone validator source is intentionally excluded from version control for now, so FASTQ reruns require the retained local companion file.

## Current validation boundary

The script and local-only validator have been syntax-checked in R, and validation-only mode has passed for both manifests, the archived processed-data contracts, and the Zenodo reference configuration. The batch-specific TSV serializers were also tested using the archived count/taxonomy matrices as stand-ins for the unavailable DADA2 objects; every generated canonical TSV/CSV matched its archived counterpart line-for-line. This tests file schemas, filtering, identifiers, ordering, taxonomy serialization, lineage formatting, and relative-abundance calculations without claiming that denoising itself was rerun.

The denoising and taxonomy run cannot be tested until the raw source sequencing data deposited in the SRA and the accession-to-sample mapping are publicly available. At that time, review quality profiles, read retention, paired-read overlap, and the automatic archived-output comparison before replacing any archived product. Change trimming values only if the released reads demonstrate that the deposited data differ from the historical inputs.

This workflow uses R, DADA2, the `biomformat` interface, the BIOM format, the SILVA database, the versioned SILVA 138.1 DADA2 training dataset, and the Sequence Read Archive [@r-core; @dada2; @biomformat; @biom; @silva; @silva-dada2-138-1; @sra].

## References

- <a id="ref-dada2"></a>Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods, 13*(7), 581–583. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- <a id="ref-sra"></a>Leinonen, R., Sugawara, H., Shumway, M., & International Nucleotide Sequence Database Collaboration. (2011). The Sequence Read Archive. *Nucleic Acids Research, 39*(Database issue), D19–D21. [https://doi.org/10.1093/nar/gkq1019](https://doi.org/10.1093/nar/gkq1019)
- <a id="ref-biom"></a>McDonald, D., Clemente, J. C., Kuczynski, J., Rideout, J. R., Stombaugh, J., Wendel, D., Wilke, A., Huse, S., Hufnagle, J., Meyer, F., Knight, R., & Caporaso, J. G. (2012). The Biological Observation Matrix (BIOM) format or: How I learned to stop worrying and love the ome-ome. *GigaScience, 1*(1), 7. [https://doi.org/10.1186/2047-217X-1-7](https://doi.org/10.1186/2047-217X-1-7)
- <a id="ref-silva-dada2-138-1"></a>McLaren, M. R., & Callahan, B. J. (2021). *Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2* (Version 1) [Data set]. Zenodo. [https://doi.org/10.5281/zenodo.4587955](https://doi.org/10.5281/zenodo.4587955)
- <a id="ref-biomformat"></a>McMurdie, P. J., & Paulson, J. N. (2026). *biomformat: An interface package for the BIOM file format* (Version 1.40.0) [R/Bioconductor package]. [https://doi.org/10.18129/B9.bioc.biomformat](https://doi.org/10.18129/B9.bioc.biomformat)
- <a id="ref-silva"></a>Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. *Nucleic Acids Research, 41*(D1), D590–D596. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- <a id="ref-r-core"></a>R Core Team. (2026). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing. [https://doi.org/10.32614/R.manuals](https://doi.org/10.32614/R.manuals)
