# Generate KO-level treatment z-scores for the KEGG pathways used in the
# Figure 7, Figure 8, and Figure S8-S16 heatmaps.
#
# This is a portable replacement for calculate_ZScore_KO2kegg() and
# write_KOs_ZScore() in scripts/legacy_scripts/picrust_legacy.R.
# It uses only canonical PICRUSt2 KO predictions, the canonical sample manifest,
# and the offline KEGG membership reference installed with ggpicrust2. Legacy CSV
# exports are not read here; they are used only by the separate validator.

# Load Packages ----
suppressPackageStartupMessages({
  library(stats)
  library(tools)
  library(utils)
})

ko_z_score_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
ko_z_score_reference_relative_path <- file.path(
  "data", "reference", "ggpicrust2", "kegg_reference.RData"
)
ko_z_score_reference_md5 <- "67728b7ea0ca954b767f77406e7a44b9"

# Explicit figure-source specification. These scientific parameters identify
# which cohort/pathway tables to generate; no checkpoint filenames are parsed.
ko_z_score_targets <- data.frame(
  batch_id = c(
    rep("batch01_iontorrent_dams-f0_f1-pnd21", 15L),
    rep("batch02_illumina-miseq_f1-f2_pnd365", 4L)
  ),
  generation = c(rep("F1", 17L), rep("F2", 2L)),
  collection_day_pnd = c(rep(21L, 15L), rep(365L, 4L)),
  cohort_label = c(rep("f1_pnd21", 15L), rep("f1_pnd365", 2L), rep("f2_pnd365", 2L)),
  pathway_id = c(
    "ko00190", "ko00310", "ko00430", "ko00480", "ko00514",
    "ko00623", "ko00740", "ko00760", "ko00930", "ko00944",
    "ko03050", "ko03450", "ko04011", "ko04962", "ko05219",
    "ko00965", "ko03050",
    "ko00311", "ko00965"
  ),
  stringsAsFactors = FALSE
)

find_ko_z_score_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata", "sample_manifest_all_batches.csv"),
      file.path(path, "results", "functional_prediction"),
      file.path(path, "scripts", "04_functional_prediction")
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the New_repo project root.", call. = FALSE)
    }
    path <- parent
  }
}

is_ko_z_score_primary_sample <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

load_kegg_ko_reference <- function(project_root) {
  reference_file <- file.path(project_root, ko_z_score_reference_relative_path)
  if (!file.exists(reference_file)) {
    stop(
      "Missing repository KEGG reference: ", reference_file,
      ". Restore data/reference/ggpicrust2/kegg_reference.RData from the repository.",
      call. = FALSE
    )
  }
  observed_md5 <- unname(tools::md5sum(reference_file))
  if (!identical(observed_md5, ko_z_score_reference_md5)) {
    stop(
      "Repository KEGG reference checksum mismatch: ", observed_md5,
      ". Expected: ", ko_z_score_reference_md5, ".",
      call. = FALSE
    )
  }

  reference_environment <- new.env(parent = emptyenv())
  loaded_names <- load(reference_file, envir = reference_environment)
  if (!identical(loaded_names, "ko_to_kegg_reference")) {
    stop("The ggpicrust2 KEGG reference has an unexpected structure.", call. = FALSE)
  }

  reference <- get("ko_to_kegg_reference", envir = reference_environment)
  reference <- as.data.frame(reference, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(dim(reference), c(306L, 326L)) || names(reference)[1L] != "X1") {
    stop("The ggpicrust2 KEGG reference is missing its pathway-ID column.", call. = FALSE)
  }
  if (anyDuplicated(reference$X1)) {
    stop("The ggpicrust2 KEGG reference contains duplicated pathway IDs.", call. = FALSE)
  }

  reference
}

get_pathway_ko_ids <- function(reference, pathway_id) {
  pathway_row <- reference[reference$X1 == pathway_id, , drop = FALSE]
  if (nrow(pathway_row) != 1L) {
    stop("Expected one KEGG reference row for pathway ", pathway_id, ".", call. = FALSE)
  }

  ko_ids <- as.character(unlist(pathway_row[1L, -1L, drop = FALSE], use.names = FALSE))
  ko_ids <- unique(ko_ids[!is.na(ko_ids) & nzchar(ko_ids)])
  if (length(ko_ids) == 0L || any(!grepl("^K[0-9]{5}$", ko_ids))) {
    stop("Invalid KEGG Ortholog membership for pathway ", pathway_id, ".", call. = FALSE)
  }

  ko_ids
}

read_ko_z_score_cohort <- function(project_root, target) {
  abundance_file <- file.path(
    project_root,
    "results", "functional_prediction", target$batch_id,
    "kegg_ortholog_abundance.tsv"
  )
  metadata_file <- file.path(
    project_root, "data", "metadata", "sample_manifest_all_batches.csv"
  )
  if (!file.exists(abundance_file)) {
    stop("Missing canonical KO abundance input: ", abundance_file, call. = FALSE)
  }
  if (!file.exists(metadata_file)) {
    stop("Missing canonical sample manifest: ", metadata_file, call. = FALSE)
  }

  abundance <- utils::read.delim(
    abundance_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  metadata <- utils::read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )

  required_metadata <- c(
    "sample_id", "batch_id", "generation", "collection_day_pnd",
    "treatment", "include_primary_analysis"
  )
  if (!all(required_metadata %in% names(metadata))) {
    stop("The canonical manifest is missing required cohort columns.", call. = FALSE)
  }
  if (ncol(abundance) < 2L || names(abundance)[1L] != "ID") {
    stop("The canonical KO table must start with an ID column.", call. = FALSE)
  }
  abundance$ID <- sub("^ko:", "", abundance$ID)

  included <-
    metadata$batch_id == target$batch_id &
    metadata$generation == target$generation &
    suppressWarnings(as.integer(metadata$collection_day_pnd)) == target$collection_day_pnd &
    is_ko_z_score_primary_sample(metadata$include_primary_analysis)
  included[is.na(included)] <- FALSE
  metadata <- metadata[included, , drop = FALSE]
  metadata$treatment <- factor(
    metadata$treatment,
    levels = ko_z_score_treatment_levels
  )

  if (nrow(metadata) == 0L || anyNA(metadata$treatment)) {
    stop("No valid samples and treatment labels remain for ", target$cohort_label, ".", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("The canonical manifest contains duplicated sample IDs.", call. = FALSE)
  }
  treatment_counts <- table(metadata$treatment)
  if (any(treatment_counts == 0L)) {
    stop("At least one treatment is absent from ", target$cohort_label, ".", call. = FALSE)
  }

  missing_samples <- setdiff(metadata$sample_id, names(abundance))
  if (length(missing_samples) > 0L) {
    stop(
      "Included samples missing from the legacy KO2 table: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(abundance$ID) || any(!grepl("^K[0-9]{5}$", abundance$ID))) {
    stop("The canonical KO table contains duplicated or malformed KO IDs.", call. = FALSE)
  }

  sample_ids <- metadata$sample_id
  numeric_values <- as.matrix(abundance[, sample_ids, drop = FALSE])
  storage.mode(numeric_values) <- "double"
  if (any(!is.finite(numeric_values)) || any(numeric_values < 0)) {
    stop("The canonical KO table contains invalid abundance values.", call. = FALSE)
  }
  abundance[, sample_ids] <- numeric_values

  list(
    abundance = abundance[, c("ID", sample_ids), drop = FALSE],
    metadata = metadata,
    sample_ids = sample_ids,
    treatment_counts = treatment_counts,
    source_files = c(abundance = abundance_file, metadata = metadata_file)
  )
}

calculate_pathway_ko_z_scores <- function(
    abundance,
    metadata,
    pathway_ko_ids,
    pathway_id) {
  selected <- abundance$ID %in% pathway_ko_ids
  pathway_abundance <- abundance[selected, , drop = FALSE]
  if (nrow(pathway_abundance) == 0L) {
    stop("No measured KO members remain for pathway ", pathway_id, ".", call. = FALSE)
  }

  sample_ids <- metadata$sample_id
  values <- as.matrix(pathway_abundance[, sample_ids, drop = FALSE])
  storage.mode(values) <- "double"
  group_means <- do.call(cbind, lapply(ko_z_score_treatment_levels, function(treatment) {
    group_values <- t(values[, metadata$treatment == treatment, drop = FALSE])
    colMeans(group_values)
  }))
  colnames(group_means) <- ko_z_score_treatment_levels

  row_mean <- apply(group_means, 1L, mean)
  row_sd <- apply(group_means, 1L, stats::sd)
  if (any(!is.finite(row_sd))) {
    bad_ids <- pathway_abundance$ID[!is.finite(row_sd)]
    stop(
      "Cannot calculate treatment variability for KO(s): ",
      paste(bad_ids, collapse = ", "),
      call. = FALSE
    )
  }

  z_scores <- sweep(group_means, 1L, row_mean, "-")
  z_scores <- sweep(z_scores, 1L, row_sd, "/")
  # The legacy export retained pathway members with identical group means as
  # explicit NA rows; downstream heatmaps omitted them with na.omit().
  z_scores[row_sd == 0, ] <- NA_real_
  output <- data.frame(
    ko_id = pathway_abundance$ID,
    z_scores,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(output) <- c("ko_id", ko_z_score_treatment_levels)
  output
}

ko_z_score_output_file <- function(output_root, target) {
  file.path(
    output_root,
    target$batch_id,
    "kegg_pathway_ko_z_scores",
    paste0(target$cohort_label, "_", target$pathway_id, "_ko_z_scores.csv")
  )
}

generate_kegg_pathway_ko_z_scores <- function(
    project_root = find_ko_z_score_project_root(),
    output_root = file.path(project_root, "results", "functional_prediction"),
    targets = ko_z_score_targets) {
  required_target_columns <- c(
    "batch_id", "generation", "collection_day_pnd", "cohort_label",
    "pathway_id"
  )
  if (!all(required_target_columns %in% names(targets)) || nrow(targets) == 0L) {
    stop("The z-score target specification is incomplete.", call. = FALSE)
  }
  target_keys <- paste(targets$cohort_label, targets$pathway_id, sep = "/")
  if (anyDuplicated(target_keys)) {
    stop("The z-score target specification contains duplicates.", call. = FALSE)
  }

  reference <- load_kegg_ko_reference(project_root)
  cohort_keys <- paste(
    targets$batch_id, targets$generation, targets$collection_day_pnd, sep = "/"
  )
  cohorts <- lapply(unique(cohort_keys), function(key) {
    target <- targets[match(key, cohort_keys), , drop = FALSE]
    read_ko_z_score_cohort(project_root, target)
  })
  names(cohorts) <- unique(cohort_keys)

  generated <- lapply(seq_len(nrow(targets)), function(index) {
    target <- targets[index, , drop = FALSE]
    cohort <- cohorts[[cohort_keys[[index]]]]
    pathway_ko_ids <- get_pathway_ko_ids(reference, target$pathway_id)
    z_scores <- calculate_pathway_ko_z_scores(
      abundance = cohort$abundance,
      metadata = cohort$metadata,
      pathway_ko_ids = pathway_ko_ids,
      pathway_id = target$pathway_id
    )

    output_file <- ko_z_score_output_file(output_root, target)
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      z_scores,
      file = output_file,
      row.names = FALSE,
      na = "NA",
      fileEncoding = "UTF-8"
    )

    data.frame(
      batch_id = target$batch_id,
      generation = target$generation,
      collection_day_pnd = target$collection_day_pnd,
      pathway_id = target$pathway_id,
      measured_ko_count = nrow(z_scores),
      output_file = normalizePath(output_file, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  })

  manifest <- do.call(rbind, generated)
  rownames(manifest) <- NULL
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  generated <- generate_kegg_pathway_ko_z_scores()
  message(
    "Generated ", nrow(generated),
    " KEGG-pathway KO z-score source tables (", sum(generated$measured_ko_count),
    " pathway-KO rows)."
  )
  print(generated, row.names = FALSE)
}
