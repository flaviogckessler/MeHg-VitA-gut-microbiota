# DADA2 preprocessing for the two 16S rRNA sequencing batches.
#
# The raw/source sequencing reads are deposited in the NCBI Sequence Read
# Archive (SRA), BioProject PRJNA1489663. They are publicly available but are
# not distributed in this repository.
#
# Run this script from the repository root with Rscript. Set:
#   DADA2_BATCH_ID   one batch identifier listed below
#   RAW_FASTQ_DIR    directory containing demultiplexed FASTQ files
# Optional:
#   SILVA_TAXONOMY_FASTA  path to silva_nr99_v138.1_train_set.fa.gz
#   SILVA_SPECIES_FASTA   path to silva_species_assignment_v138.1.fa.gz
#   DADA2_OUTPUT_DIR    rerun output directory
#   DADA2_VALIDATE_ONLY true to validate configuration and archived products

# Load Packages ----
preprocess_validate_only <-
  tolower(Sys.getenv("DADA2_VALIDATE_ONLY", unset = "false")) == "true"
suppressPackageStartupMessages(library(tools))
if (!preprocess_validate_only) {
  suppressPackageStartupMessages({
    library(dada2)
    library(biomformat)
  })
}

sra_bioproject <- "PRJNA1489663"
silva_zenodo_doi <- "10.5281/zenodo.4587955"
silva_taxonomy_md5 <- "6b41db7139834c71171f8ce5b5918fc6"
silva_species_md5 <- "f21c2d97c79ff07c17949a9622371a4c"
batch01 <- "batch01_iontorrent_dams-f0_f1-pnd21"
batch02 <- "batch02_illumina-miseq_f1-f2_pnd365"
batch_id <- Sys.getenv("DADA2_BATCH_ID", unset = "")
validate_only <- preprocess_validate_only

script_option <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_option) == 1L) {
  normalizePath(
    file.path(dirname(sub("^--file=", "", script_option)), "..", ".."),
    mustWork = TRUE
  )
} else {
  normalizePath(".", mustWork = TRUE)
}

if (!batch_id %in% c(batch01, batch02)) {
  stop(
    "Set DADA2_BATCH_ID to either ", batch01, " or ", batch02, ".",
    call. = FALSE
  )
}

manifest_file <- file.path(repo_root, "data", "metadata", batch_id, "sample_manifest.csv")
canonical_dir <- file.path(repo_root, "data", "processed", batch_id)
validator_file <- file.path(
  repo_root, "scripts", "01_sequence_processing", "validate_preprocess_16s_output.R"
)
if (!file.exists(validator_file)) {
  stop("Preprocessing validator not found: ", validator_file, call. = FALSE)
}
source(validator_file, local = TRUE)
silva_taxonomy_file <- Sys.getenv(
  "SILVA_TAXONOMY_FASTA",
  unset = file.path(
    repo_root, "data", "external", "silva", "silva_nr99_v138.1_train_set.fa.gz"
  )
)
silva_species_file <- Sys.getenv(
  "SILVA_SPECIES_FASTA",
  unset = file.path(
    repo_root, "data", "external", "silva", "silva_species_assignment_v138.1.fa.gz"
  )
)

if (!file.exists(manifest_file)) {
  stop("Canonical sample manifest not found: ", manifest_file, call. = FALSE)
}

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
if (!"sample_id" %in% names(manifest) || anyDuplicated(manifest$sample_id)) {
  stop("The canonical manifest must contain unique sample_id values.", call. = FALSE)
}

if (validate_only) {
  validate_archived_preprocess_contract(batch_id, canonical_dir, manifest)
  reference_files <- c(silva_taxonomy_file, silva_species_file)
  if (all(file.exists(reference_files))) {
    observed_md5 <- tolower(unname(tools::md5sum(reference_files)))
    expected_md5 <- c(silva_taxonomy_md5, silva_species_md5)
    if (!identical(observed_md5, expected_md5)) {
      stop(
        "A local SILVA FASTA checksum does not match Zenodo record 4587955.",
        call. = FALSE
      )
    }
    reference_status <- "the two local SILVA FASTA checksums match Zenodo"
  } else {
    reference_status <- paste0(
      "the Zenodo DOI/checksum configuration is present; local FASTAs are not yet downloaded"
    )
  }
  message(
    "Validation passed for ", batch_id, ": the script parsed, the canonical ",
    "manifest and archived processed-data contract are consistent, and ",
    reference_status, "."
  )
  quit(save = "no", status = 0L)
}

reference_files <- c(silva_taxonomy_file, silva_species_file)
if (!all(file.exists(reference_files))) {
  stop(
    "SILVA 138.1 DADA2 reference FASTAs are missing. Download the taxonomy and ",
    "species files from https://doi.org/", silva_zenodo_doi,
    " as described in data/external/README.md.",
    call. = FALSE
  )
}
observed_md5 <- tolower(unname(tools::md5sum(reference_files)))
expected_md5 <- c(silva_taxonomy_md5, silva_species_md5)
if (!identical(observed_md5, expected_md5)) {
  stop(
    "A local SILVA FASTA checksum does not match Zenodo record 4587955.",
    call. = FALSE
  )
}

raw_fastq_dir <- Sys.getenv("RAW_FASTQ_DIR", unset = "")
if (!nzchar(raw_fastq_dir) || !dir.exists(raw_fastq_dir)) {
  stop(
    "Raw FASTQ files are not distributed in this repository. The source ",
    "sequencing reads are deposited in NCBI SRA under BioProject ",
    sra_bioproject, ". Download and demultiplex the public reads, verify their ",
    "sample assignments, and set RAW_FASTQ_DIR to the resulting directory.",
    call. = FALSE
  )
}
raw_fastq_dir <- normalizePath(raw_fastq_dir, mustWork = TRUE)

required_packages <- c("dada2", "biomformat")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the required Bioconductor packages before processing: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

output_dir <- Sys.getenv(
  "DADA2_OUTPUT_DIR",
  unset = file.path(repo_root, "data", "processed", batch_id, "dada2_rerun")
)
if (normalizePath(output_dir, mustWork = FALSE) == normalizePath(canonical_dir, mustWork = TRUE)) {
  stop(
    "DADA2_OUTPUT_DIR cannot be the archived processed-data directory. ",
    "Write to a review directory and validate the rerun first.",
    call. = FALSE
  )
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
filtered_dir <- file.path(output_dir, "filtered_fastq")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

if (batch_id == batch01) {
  fastq <- sort(list.files(
    raw_fastq_dir,
    pattern = "[.](fastq|fq)([.]gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  ))
  if (!length(fastq)) {
    stop("No single-end FASTQ files were found in RAW_FASTQ_DIR.", call. = FALSE)
  }

  sample_ids <- sub(
    "[.](fastq|fq)([.]gz)?$", "", basename(fastq), ignore.case = TRUE
  )
  sample_ids <- sub("[.]IonXpress.*$", "", sample_ids)
  if (anyDuplicated(sample_ids)) {
    stop("FASTQ filenames do not produce unique sample identifiers.", call. = FALSE)
  }
  if (!setequal(sample_ids, manifest$sample_id)) {
    stop(
      "FASTQ sample identifiers must match the canonical manifest. Missing: ",
      paste(setdiff(manifest$sample_id, sample_ids), collapse = ", "),
      "; unexpected: ",
      paste(setdiff(sample_ids, manifest$sample_id), collapse = ", "),
      ". See scripts/01_sequence_processing/README.md.",
      call. = FALSE
    )
  }

  filtered <- file.path(filtered_dir, paste0(sample_ids, "_filtered.fastq.gz"))
  names(filtered) <- sample_ids

  filtering <- dada2::filterAndTrim(
    fastq, filtered,
    truncLen = 275, trimLeft = 20,
    maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
    compress = TRUE, multithread = TRUE, verbose = TRUE
  )
  errors <- dada2::learnErrors(filtered, multithread = TRUE)
  denoised <- dada2::dada(
    filtered, err = errors, pool = "pseudo",
    HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32,
    multithread = TRUE
  )
  sequence_table <- dada2::makeSequenceTable(denoised)
  sequence_table <- dada2::removeBimeraDenovo(
    sequence_table, method = "consensus", multithread = TRUE, verbose = TRUE
  )

  read_tracking <- data.frame(
    sample_id = sample_ids,
    input = filtering[, "reads.in"],
    filtered = filtering[, "reads.out"],
    denoised = vapply(
      denoised, function(x) sum(dada2::getUniques(x)), numeric(1)
    ),
    nonchim = rowSums(sequence_table),
    check.names = FALSE
  )
  parameter_values <- c(
    platform = "Ion Torrent; single-end",
    truncLen = "275",
    trimLeft = "20",
    maxEE = "2",
    pooling = "pseudo",
    homopolymer_gap_penalty = "-1",
    band_size = "32"
  )
} else {
  forward <- sort(list.files(
    raw_fastq_dir,
    pattern = "_R1_001[.](fastq|fq)([.]gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  ))
  reverse <- sort(list.files(
    raw_fastq_dir,
    pattern = "_R2_001[.](fastq|fq)([.]gz)?$",
    full.names = TRUE,
    ignore.case = TRUE
  ))
  if (!length(forward) || length(forward) != length(reverse)) {
    stop("Paired FASTQ files are absent or forward/reverse counts differ.", call. = FALSE)
  }

  forward_ids <- sub(
    "_R1_001[.](fastq|fq)([.]gz)?$", "", basename(forward), ignore.case = TRUE
  )
  reverse_ids <- sub(
    "_R2_001[.](fastq|fq)([.]gz)?$", "", basename(reverse), ignore.case = TRUE
  )
  if (!identical(forward_ids, reverse_ids) || anyDuplicated(forward_ids)) {
    stop("Forward/reverse FASTQ filenames are not uniquely paired.", call. = FALSE)
  }
  sample_ids <- forward_ids
  expected_sample_ids <- c(manifest$sample_id, "t3")
  if (!setequal(sample_ids, expected_sample_ids)) {
    stop(
      "To reproduce the archived batch 02 products, FASTQ identifiers must ",
      "match the canonical manifest plus historical sample t3. Missing: ",
      paste(setdiff(expected_sample_ids, sample_ids), collapse = ", "),
      "; unexpected: ",
      paste(setdiff(sample_ids, expected_sample_ids), collapse = ", "),
      ". See scripts/01_sequence_processing/README.md.",
      call. = FALSE
    )
  }

  filtered_forward <- file.path(
    filtered_dir, paste0(sample_ids, "_F_filtered.fastq.gz")
  )
  filtered_reverse <- file.path(
    filtered_dir, paste0(sample_ids, "_R_filtered.fastq.gz")
  )
  names(filtered_forward) <- sample_ids
  names(filtered_reverse) <- sample_ids

  filtering <- dada2::filterAndTrim(
    forward, filtered_forward, reverse, filtered_reverse,
    truncLen = c(149, 70), trimLeft = c(20, 0),
    maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
    compress = TRUE, matchIDs = TRUE, multithread = TRUE, verbose = TRUE
  )
  errors_forward <- dada2::learnErrors(filtered_forward, multithread = TRUE)
  errors_reverse <- dada2::learnErrors(filtered_reverse, multithread = TRUE)
  denoised_forward <- dada2::dada(
    filtered_forward, err = errors_forward, pool = "pseudo", multithread = TRUE
  )
  denoised_reverse <- dada2::dada(
    filtered_reverse, err = errors_reverse, pool = "pseudo", multithread = TRUE
  )
  merged <- dada2::mergePairs(
    denoised_forward, filtered_forward,
    denoised_reverse, filtered_reverse,
    verbose = TRUE
  )
  sequence_table <- dada2::makeSequenceTable(merged)
  sequence_table <- dada2::removeBimeraDenovo(
    sequence_table, method = "consensus", multithread = TRUE, verbose = TRUE
  )

  read_tracking <- data.frame(
    sample_id = sample_ids,
    input = filtering[, "reads.in"],
    filtered = filtering[, "reads.out"],
    denoised_forward = vapply(
      denoised_forward, function(x) sum(dada2::getUniques(x)), numeric(1)
    ),
    denoised_reverse = vapply(
      denoised_reverse, function(x) sum(dada2::getUniques(x)), numeric(1)
    ),
    merged = vapply(merged, function(x) sum(dada2::getUniques(x)), numeric(1)),
    nonchim = rowSums(sequence_table),
    check.names = FALSE
  )
  parameter_values <- c(
    platform = "Illumina MiSeq; paired-end",
    truncLen = "149,70",
    trimLeft = "20,0",
    maxEE = "2,2",
    pooling = "pseudo",
    paired_id_check = "true",
    archived_sample_compatibility = "canonical manifest plus t3",
    relative_taxa_filter = "Bacteria; Phylum/Class/Family assigned; prevalence > 6; total > 2; maximum > 1",
    relative_id_scheme = "SV%05d based on unfiltered ASV position"
  )
}

if (!ncol(sequence_table)) {
  stop("No non-chimeric ASVs remained after preprocessing.", call. = FALSE)
}
if (any(rowSums(sequence_table) == 0)) {
  stop("At least one sample has zero non-chimeric reads.", call. = FALSE)
}
rownames(sequence_table) <- sample_ids

taxonomy <- dada2::assignTaxonomy(
  sequence_table,
  refFasta = silva_taxonomy_file,
  multithread = TRUE,
  verbose = TRUE
)
taxonomy <- dada2::addSpecies(
  taxonomy,
  refFasta = silva_species_file,
  verbose = TRUE
)

if (batch_id == batch01) {
  counts_table <- data.frame(
    OTU = colnames(sequence_table), t(sequence_table), check.names = FALSE
  )
  taxonomy_table <- data.frame(
    OTU = rownames(taxonomy), taxonomy, row.names = NULL, check.names = FALSE
  )
  lineage <- apply(
    taxonomy, 1L,
    function(x) paste0(
      c("k__", "p__", "c__", "o__", "f__", "g__", "s__"),
      ifelse(is.na(x), "NA", x), collapse = ";"
    )
  )
  counts_with_taxonomy <- cbind(
    taxonomy_table,
    lineage = lineage,
    counts_table[, -1, drop = FALSE]
  )
  counts_with_taxonomy$lineage <- paste0('"', counts_with_taxonomy$lineage, '"')
  renamed_taxonomy <- data.frame(
    ASV = paste0("ASV", seq_len(nrow(taxonomy))),
    taxonomy, row.names = NULL, check.names = FALSE
  )

  write.table(
    counts_table, file.path(output_dir, "asv_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(
    taxonomy_table, file.path(output_dir, "asv_taxonomy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(
    counts_with_taxonomy,
    file.path(output_dir, "asv_counts_with_taxonomy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.table(
    renamed_taxonomy,
    file.path(output_dir, "asv_taxonomy_renamed_ids.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "NA"
  )
  tracking_table <- read_tracking[, -1, drop = FALSE]
  rownames(tracking_table) <- read_tracking$sample_id
  write.csv(
    tracking_table,
    file.path(output_dir, "read_tracking.csv"),
    row.names = TRUE, quote = TRUE
  )
} else {
  write.table(
    t(sequence_table), file.path(output_dir, "asv_counts.tsv"),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )
  taxonomy_text <- apply(
    taxonomy, 1L,
    function(x) paste(ifelse(is.na(x), "NA", x), collapse = ";")
  )
  write.table(
    data.frame(sequence = rownames(taxonomy), taxonomy = taxonomy_text),
    file.path(output_dir, "asv_taxonomy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )

  # Retain mitochondrial-family assignments intentionally. The archived
  # batch-02 tables used downstream include them, and the historical plural
  # predicate (Family != "Mitochondrias") never excluded the actual singular
  # Family value "Mitochondria".
  keep <- !is.na(taxonomy[, "Kingdom"]) & taxonomy[, "Kingdom"] == "Bacteria" &
    !is.na(taxonomy[, "Phylum"]) & !is.na(taxonomy[, "Class"]) &
    !is.na(taxonomy[, "Family"]) &
    colSums(sequence_table > 0) > 6 & colSums(sequence_table) > 2 &
    apply(sequence_table, 2L, max) > 1
  filtered_table <- sequence_table[, keep, drop = FALSE]
  relative_table <- sweep(filtered_table, 1L, rowSums(filtered_table), "/")
  relative_taxonomy <- taxonomy[keep, 1:6, drop = FALSE]
  # Preserve the label used by the archived filtered table without altering
  # the canonical full taxonomy export.
  legacy_genus <- !is.na(relative_taxonomy[, "Genus"]) &
    relative_taxonomy[, "Family"] == "Ruminococcaceae" &
    relative_taxonomy[, "Genus"] == "Incertae Sedis"
  relative_taxonomy[legacy_genus, "Genus"] <- "Ruminococcaceae Incertae Sedis"
  relative_with_taxonomy <- data.frame(
    ASV = sprintf("SV%05d", which(keep)),
    relative_taxonomy,
    t(relative_table),
    row.names = NULL,
    check.names = FALSE
  )
  write.table(
    relative_with_taxonomy,
    file.path(output_dir, "asv_relative_abundance_with_taxonomy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
  write.csv(
    read_tracking,
    file.path(output_dir, "read_tracking.csv"),
    row.names = FALSE, quote = FALSE
  )
}

dada2::uniquesToFasta(
  sequence_table,
  fout = file.path(output_dir, "representative_sequences.fasta"),
  ids = colnames(sequence_table)
)
feature_table <- biomformat::make_biom(data = t(sequence_table))
biomformat::write_biom(
  feature_table,
  biom_file = file.path(output_dir, "feature_table.biom")
)

parameter_table <- data.frame(
  parameter = c(
    "batch_id", "sra_bioproject", "silva_zenodo_doi",
    "silva_taxonomy_file", "silva_taxonomy_md5",
    "silva_species_file", "silva_species_md5",
    names(parameter_values)
  ),
  value = c(
    batch_id, sra_bioproject, silva_zenodo_doi,
    normalizePath(silva_taxonomy_file), silva_taxonomy_md5,
    normalizePath(silva_species_file), silva_species_md5,
    parameter_values
  ),
  row.names = NULL
)
write.table(
  parameter_table,
  file.path(output_dir, "preprocessing_parameters.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "preprocessing_session_info.txt")
)

validate_preprocess_16s_output(batch_id, output_dir, canonical_dir)
message("All generated canonical products match the archived processed data.")
message("Preprocessing completed. Review outputs in: ", normalizePath(output_dir))
