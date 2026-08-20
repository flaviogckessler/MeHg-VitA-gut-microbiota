# Validation helpers for the archived and newly generated DADA2 products.

# Load Packages ----
suppressPackageStartupMessages(library(utils))
if (tolower(Sys.getenv("DADA2_VALIDATE_ONLY", unset = "false")) != "true") {
  suppressPackageStartupMessages(library(biomformat))
}

same_values <- function(x, y, tolerance = 0) {
  if (!identical(dim(x), dim(y)) || !identical(dimnames(x), dimnames(y))) {
    return(FALSE)
  }
  if (is.numeric(x) && is.numeric(y)) {
    return(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)))
  }
  identical(is.na(x), is.na(y)) && all(x[!is.na(x)] == y[!is.na(y)])
}

read_sequence_fasta <- function(path) {
  fasta <- readLines(path, warn = FALSE)
  if (length(fasta) %% 2L != 0L || any(!startsWith(fasta[seq(1L, length(fasta), 2L)], ">"))) {
    stop("Expected a two-line-per-record FASTA file: ", path, call. = FALSE)
  }
  list(
    id = substring(fasta[seq(1L, length(fasta), 2L)], 2L),
    sequence = fasta[seq(2L, length(fasta), 2L)]
  )
}

read_asv_counts <- function(path) {
  as.matrix(read.delim(
    path,
    row.names = 1L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
}

validate_archived_preprocess_contract <- function(batch_id, canonical_dir, manifest) {
  batch01 <- "batch01_iontorrent_dams-f0_f1-pnd21"
  batch02 <- "batch02_illumina-miseq_f1-f2_pnd365"
  required <- if (batch_id == batch01) {
    c(
      "asv_counts.tsv", "asv_taxonomy.tsv", "asv_counts_with_taxonomy.tsv",
      "asv_taxonomy_renamed_ids.tsv", "representative_sequences.fasta",
      "feature_table.biom", "read_tracking.csv"
    )
  } else {
    c(
      "asv_counts.tsv", "asv_taxonomy.tsv",
      "asv_relative_abundance_with_taxonomy.tsv",
      "representative_sequences.fasta", "feature_table.biom"
    )
  }
  missing <- required[!file.exists(file.path(canonical_dir, required))]
  if (length(missing)) {
    stop("Archived preprocessing product(s) missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  counts <- read_asv_counts(file.path(canonical_dir, "asv_counts.tsv"))
  storage.mode(counts) <- "numeric"
  sequences <- rownames(counts)
  fasta <- read_sequence_fasta(file.path(canonical_dir, "representative_sequences.fasta"))
  if (!identical(fasta$id, sequences) || !identical(fasta$sequence, sequences)) {
    stop("Archived representative FASTA does not match the ASV count-table order.", call. = FALSE)
  }

  if (batch_id == batch01) {
    if (nrow(counts) != 1701L || ncol(counts) != 54L || any(nchar(sequences) != 255L)) {
      stop("Batch 01 archived count-table dimensions or sequence lengths changed.", call. = FALSE)
    }
    if (!setequal(colnames(counts), manifest$sample_id)) {
      stop("Batch 01 archived samples do not match the canonical manifest.", call. = FALSE)
    }

    taxonomy <- read.delim(
      file.path(canonical_dir, "asv_taxonomy.tsv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    combined <- read.delim(
      file.path(canonical_dir, "asv_counts_with_taxonomy.tsv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    renamed <- read.delim(
      file.path(canonical_dir, "asv_taxonomy_renamed_ids.tsv"),
      header = FALSE, check.names = FALSE, stringsAsFactors = FALSE
    )
    tracking <- read.csv(
      file.path(canonical_dir, "read_tracking.csv"),
      row.names = 1L, check.names = FALSE
    )
    expected_lineage <- apply(
      taxonomy[, 2:8, drop = FALSE], 1L,
      function(x) paste0(
        c("k__", "p__", "c__", "o__", "f__", "g__", "s__"),
        ifelse(is.na(x), "NA", x), collapse = ";"
      )
    )
    if (!identical(taxonomy[[1L]], sequences) ||
        !same_values(as.matrix(combined[, 1:8]), as.matrix(taxonomy)) ||
        !identical(combined$lineage, expected_lineage) ||
        !identical(colnames(combined)[-(1:9)], colnames(counts)) ||
        !same_values(unname(as.matrix(combined[, -(1:9)])), unname(counts)) ||
        !identical(renamed[[1L]], paste0("ASV", seq_len(nrow(counts)))) ||
        !same_values(unname(as.matrix(renamed[, -1])), unname(as.matrix(taxonomy[, -1]))) ||
        !identical(rownames(tracking), colnames(counts))) {
      stop("Batch 01 archived tables are not mutually consistent.", call. = FALSE)
    }
  } else if (batch_id == batch02) {
    if (nrow(counts) != 11704L || ncol(counts) != 49L || any(nchar(sequences) != 150L)) {
      stop("Batch 02 archived count-table dimensions or sequence lengths changed.", call. = FALSE)
    }
    if (!setequal(colnames(counts), c(manifest$sample_id, "t3"))) {
      stop("Batch 02 archived samples must equal the manifest plus historical sample t3.", call. = FALSE)
    }

    taxonomy_file <- read.delim(
      file.path(canonical_dir, "asv_taxonomy.tsv"),
      header = FALSE, check.names = FALSE, stringsAsFactors = FALSE,
      na.strings = character()
    )
    if (!identical(taxonomy_file[[1L]], sequences)) {
      stop("Batch 02 taxonomy order does not match the count table.", call. = FALSE)
    }
    taxonomy <- do.call(rbind, strsplit(taxonomy_file[[2L]], ";", fixed = TRUE))
    colnames(taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    taxonomy[taxonomy == "NA"] <- NA_character_

    keep <- !is.na(taxonomy[, "Kingdom"]) & taxonomy[, "Kingdom"] == "Bacteria" &
      !is.na(taxonomy[, "Phylum"]) & !is.na(taxonomy[, "Class"]) &
      !is.na(taxonomy[, "Family"]) & taxonomy[, "Family"] != "Mitochondrias" &
      rowSums(counts > 0) > 6 & rowSums(counts) > 2 & apply(counts, 1L, max) > 1
    relative <- read.delim(
      file.path(canonical_dir, "asv_relative_abundance_with_taxonomy.tsv"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    expected_ids <- sprintf("SV%05d", which(keep))
    expected_taxonomy <- taxonomy[keep, 1:6, drop = FALSE]
    legacy_genus <- !is.na(expected_taxonomy[, "Genus"]) &
      expected_taxonomy[, "Family"] == "Ruminococcaceae" &
      expected_taxonomy[, "Genus"] == "Incertae Sedis"
    expected_taxonomy[legacy_genus, "Genus"] <- "Ruminococcaceae Incertae Sedis"
    expected_relative <- sweep(counts[keep, , drop = FALSE], 2L, colSums(counts[keep, , drop = FALSE]), "/")
    observed_relative <- as.matrix(relative[, -(1:7), drop = FALSE])
    storage.mode(observed_relative) <- "numeric"

    if (sum(keep) != 2739L || !identical(relative[[1L]], expected_ids) ||
        !same_values(unname(as.matrix(relative[, 2:7])), unname(expected_taxonomy)) ||
        !identical(colnames(observed_relative), colnames(counts)) ||
        !same_values(unname(observed_relative), unname(expected_relative), tolerance = 1e-14)) {
      stop("Batch 02 archived filtered relative-abundance table is inconsistent.", call. = FALSE)
    }
  } else {
    stop("Unknown preprocessing batch: ", batch_id, call. = FALSE)
  }

  invisible(TRUE)
}

validate_preprocess_16s_output <- function(batch_id, candidate_dir, canonical_dir) {
  batch01 <- "batch01_iontorrent_dams-f0_f1-pnd21"
  expected <- if (batch_id == batch01) {
    c(
      "asv_counts.tsv", "asv_taxonomy.tsv", "asv_counts_with_taxonomy.tsv",
      "asv_taxonomy_renamed_ids.tsv", "representative_sequences.fasta",
      "feature_table.biom", "read_tracking.csv"
    )
  } else {
    c(
      "asv_counts.tsv", "asv_taxonomy.tsv",
      "asv_relative_abundance_with_taxonomy.tsv",
      "representative_sequences.fasta", "feature_table.biom"
    )
  }
  report <- data.frame(
    artifact = expected,
    status = "not_checked",
    detail = "",
    stringsAsFactors = FALSE
  )
  forbidden <- if (batch_id == batch01) {
    "asv_relative_abundance_with_taxonomy.tsv"
  } else {
    c("asv_counts_with_taxonomy.tsv", "asv_taxonomy_renamed_ids.tsv")
  }
  missing_candidate <- !file.exists(file.path(candidate_dir, expected))
  missing_canonical <- !file.exists(file.path(canonical_dir, expected))
  report$status[missing_candidate | missing_canonical] <- "FAIL"
  report$detail[missing_candidate] <- "new output is missing"
  report$detail[missing_canonical] <- "archived comparator is missing"
  unexpected <- forbidden[file.exists(file.path(candidate_dir, forbidden))]
  if (length(unexpected)) {
    report <- rbind(
      report,
      data.frame(
        artifact = unexpected,
        status = "FAIL",
        detail = "stale output from the other batch contract",
        stringsAsFactors = FALSE
      )
    )
  }

  comparable <- expected[!(missing_candidate | missing_canonical)]
  text_files <- setdiff(comparable, c("feature_table.biom", "asv_relative_abundance_with_taxonomy.tsv"))
  for (artifact in text_files) {
    i <- match(artifact, report$artifact)
    equal <- identical(
      readLines(file.path(candidate_dir, artifact), warn = FALSE),
      readLines(file.path(canonical_dir, artifact), warn = FALSE)
    )
    report$status[i] <- if (equal) "PASS" else "FAIL"
    report$detail[i] <- if (equal) "line-for-line match" else "content differs"
  }

  if ("asv_relative_abundance_with_taxonomy.tsv" %in% comparable) {
    artifact <- "asv_relative_abundance_with_taxonomy.tsv"
    i <- match(artifact, report$artifact)
    candidate <- read.delim(file.path(candidate_dir, artifact), check.names = FALSE)
    canonical <- read.delim(file.path(canonical_dir, artifact), check.names = FALSE)
    fixed_equal <- same_values(as.matrix(candidate[, 1:7]), as.matrix(canonical[, 1:7]))
    candidate_values <- as.matrix(candidate[, -(1:7), drop = FALSE])
    canonical_values <- as.matrix(canonical[, -(1:7), drop = FALSE])
    storage.mode(candidate_values) <- "numeric"
    storage.mode(canonical_values) <- "numeric"
    numeric_equal <- same_values(candidate_values, canonical_values, tolerance = 1e-14)
    equal <- fixed_equal && numeric_equal
    report$status[i] <- if (equal) "PASS" else "FAIL"
    report$detail[i] <- if (equal) "identifiers/taxonomy exact; abundance tolerance 1e-14" else "content differs"
  }

  if ("feature_table.biom" %in% comparable) {
    artifact <- "feature_table.biom"
    i <- match(artifact, report$artifact)
    if (!requireNamespace("biomformat", quietly = TRUE)) {
      report$status[i] <- "FAIL"
      report$detail[i] <- "biomformat is required for semantic BIOM comparison"
    } else {
      candidate <- as.matrix(biomformat::biom_data(
        biomformat::read_biom(file.path(candidate_dir, artifact))
      ))
      canonical <- as.matrix(biomformat::biom_data(
        biomformat::read_biom(file.path(canonical_dir, artifact))
      ))
      equal <- same_values(candidate, canonical)
      report$status[i] <- if (equal) "PASS" else "FAIL"
      report$detail[i] <- if (equal) "semantic matrix match" else "feature/sample matrix differs"
    }
  }

  report_file <- file.path(candidate_dir, "archived_output_comparison.tsv")
  write.table(report, report_file, sep = "\t", quote = FALSE, row.names = FALSE)
  if (any(report$status != "PASS")) {
    failed <- report$artifact[report$status != "PASS"]
    stop(
      "New preprocessing products do not match the archived data: ",
      paste(failed, collapse = ", "), ". See ", report_file, ".",
      call. = FALSE
    )
  }
  invisible(report)
}
