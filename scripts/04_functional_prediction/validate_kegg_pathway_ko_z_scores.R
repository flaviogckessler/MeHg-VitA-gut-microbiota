# Validate the portable KEGG-pathway KO z-score source-data workflow.
#
# Checks include regenerated-versus-committed equality, row-wise z-score
# invariants, and comparison with the unchanged CSVs exported by the legacy
# calculate_ZScore_KO2kegg() workflow.

# Load Packages ----
suppressPackageStartupMessages({
  library(stats)
  library(utils)
})

find_ko_z_score_script_directory <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidates <- c(path, file.path(path, "scripts", "04_functional_prediction"))
    for (candidate in candidates) {
      generator <- file.path(candidate, "generate_kegg_pathway_ko_z_scores.R")
      if (file.exists(generator)) return(candidate)
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_kegg_pathway_ko_z_scores.R.", call. = FALSE)
    }
    path <- parent
  }
}

script_directory <- find_ko_z_score_script_directory()
source(file.path(script_directory, "generate_kegg_pathway_ko_z_scores.R"))

check_ko_z_score <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_clean_z_score_csv <- function(file) {
  table <- utils::read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  check_ko_z_score(
    identical(names(table), c("ko_id", ko_z_score_treatment_levels)),
    paste("Unexpected columns in", file)
  )
  table
}

read_legacy_z_score_csv <- function(file) {
  table <- utils::read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  expected_legacy_columns <- c("ID", ko_z_score_treatment_levels)
  if (ncol(table) == length(expected_legacy_columns) + 1L) {
    row_index <- suppressWarnings(as.integer(table[[1L]]))
    if (identical(row_index, seq_len(nrow(table)))) {
      table <- table[, -1L, drop = FALSE]
    }
  }
  check_ko_z_score(
    identical(names(table), expected_legacy_columns),
    paste("Unexpected legacy columns in", file)
  )
  names(table)[1L] <- "ko_id"
  table
}

compare_z_score_tables <- function(observed, expected, label, tolerance = 1e-12) {
  check_ko_z_score(identical(names(observed), names(expected)), paste(label, "columns differ."))
  check_ko_z_score(identical(observed$ko_id, expected$ko_id), paste(label, "KO order differs."))

  observed_values <- as.matrix(observed[, ko_z_score_treatment_levels, drop = FALSE])
  expected_values <- as.matrix(expected[, ko_z_score_treatment_levels, drop = FALSE])
  storage.mode(observed_values) <- "double"
  storage.mode(expected_values) <- "double"
  check_ko_z_score(identical(dim(observed_values), dim(expected_values)), paste(label, "dimensions differ."))
  check_ko_z_score(
    identical(is.na(observed_values), is.na(expected_values)),
    paste(label, "NA patterns differ.")
  )

  absolute_differences <- abs(observed_values - expected_values)
  max_absolute_difference <- if (all(is.na(absolute_differences))) {
    0
  } else {
    max(absolute_differences, na.rm = TRUE)
  }
  check_ko_z_score(
    is.finite(max_absolute_difference) && max_absolute_difference <= tolerance,
    paste(label, "maximum absolute difference is", format(max_absolute_difference, scientific = TRUE))
  )
  max_absolute_difference
}

legacy_checkpoint_filename <- function(target) {
  suffixes <- c(
    f1_pnd21 = "F1_PND21",
    f1_pnd365 = "F1_1year",
    f2_pnd365 = "F2_1year"
  )
  suffix <- unname(suffixes[target$cohort_label])
  check_ko_z_score(
    length(suffix) == 1L && !is.na(suffix),
    paste("No legacy checkpoint suffix is defined for", target$cohort_label)
  )
  paste0(target$pathway_id, "__", suffix, ".csv")
}

validate_canonical_ko_inputs <- function(project_root) {
  batches <- unique(ko_z_score_targets$batch_id)
  comparisons <- lapply(batches, function(batch_id) {
    input_directory <- file.path(
      project_root, "results", "functional_prediction", batch_id
    )
    canonical <- utils::read.delim(
      file.path(input_directory, "kegg_ortholog_abundance.tsv"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    legacy <- utils::read.delim(
      file.path(input_directory, "kegg_ortholog_abundance_legacy_ko2.tsv"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    canonical$ID <- sub("^ko:", "", canonical$ID)

    check_ko_z_score(identical(dim(canonical), dim(legacy)), paste(batch_id, "KO dimensions differ."))
    check_ko_z_score(identical(names(canonical), names(legacy)), paste(batch_id, "KO columns differ."))
    check_ko_z_score(identical(canonical$ID, legacy$ID), paste(batch_id, "KO identifiers differ."))

    canonical_values <- as.matrix(canonical[, -1L, drop = FALSE])
    legacy_values <- as.matrix(legacy[, -1L, drop = FALSE])
    storage.mode(canonical_values) <- "double"
    storage.mode(legacy_values) <- "double"
    max_absolute_difference <- max(abs(canonical_values - legacy_values))
    check_ko_z_score(
      identical(max_absolute_difference, 0),
      paste(batch_id, "canonical and legacy KO values differ.")
    )

    data.frame(
      batch_id = batch_id,
      ko_rows = nrow(canonical),
      sample_columns = ncol(canonical) - 1L,
      max_absolute_difference = max_absolute_difference,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, comparisons)
}

validate_kegg_pathway_ko_z_scores <- function(
    project_root = find_ko_z_score_project_root(),
    tolerance = 1e-12) {
  generator_file <- file.path(
    project_root, "scripts", "04_functional_prediction",
    "generate_kegg_pathway_ko_z_scores.R"
  )
  generator_text <- paste(readLines(generator_file, warn = FALSE), collapse = "\n")
  forbidden_generator_dependencies <- c(
    "legacy_filename",
    "legacy_summary_tables",
    "kegg_ortholog_abundance_legacy_ko2",
    "results/functional_prediction/ZScores"
  )
  check_ko_z_score(
    !any(vapply(
      forbidden_generator_dependencies,
      grepl,
      logical(1),
      x = generator_text,
      fixed = TRUE
    )),
    "The generator contains a dependency on legacy checkpoint data."
  )
  canonical_input_comparison <- validate_canonical_ko_inputs(project_root)

  temporary_output_root <- tempfile("ko_pathway_z_score_validation_")
  dir.create(temporary_output_root, recursive = TRUE, showWarnings = FALSE)
  regenerated <- generate_kegg_pathway_ko_z_scores(
    project_root = project_root,
    output_root = temporary_output_root
  )

  legacy_candidates <- c(
    file.path(project_root, "results", "legacy_summary_tables", "figure_7_8_s8_s16"),
    file.path(project_root, "results", "functional_prediction", "ZScores")
  )
  legacy_directory <- legacy_candidates[dir.exists(legacy_candidates)][1L]
  check_ko_z_score(
    length(legacy_directory) == 1L && !is.na(legacy_directory),
    "Could not locate the legacy KO z-score checkpoint directory."
  )

  expected_output_files <- vapply(seq_len(nrow(ko_z_score_targets)), function(index) {
    ko_z_score_output_file(
      file.path(project_root, "results", "functional_prediction"),
      ko_z_score_targets[index, , drop = FALSE]
    )
  }, character(1))
  check_ko_z_score(all(file.exists(expected_output_files)), "One or more committed source tables are missing.")

  output_directories <- unique(dirname(expected_output_files))
  actual_output_files <- sort(unlist(lapply(output_directories, function(directory) {
    list.files(directory, pattern = "_ko_z_scores[.]csv$", full.names = TRUE)
  }), use.names = FALSE))
  check_ko_z_score(
    setequal(normalizePath(actual_output_files, winslash = "/"),
             normalizePath(expected_output_files, winslash = "/")),
    "The generated-source directories contain missing or unexpected CSV files."
  )

  comparisons <- lapply(seq_len(nrow(ko_z_score_targets)), function(index) {
    target <- ko_z_score_targets[index, , drop = FALSE]
    regenerated_file <- regenerated$output_file[[index]]
    committed_file <- expected_output_files[[index]]
    legacy_file <- file.path(legacy_directory, legacy_checkpoint_filename(target))
    check_ko_z_score(file.exists(legacy_file), paste("Missing legacy checkpoint:", legacy_file))

    regenerated_table <- read_clean_z_score_csv(regenerated_file)
    committed_table <- read_clean_z_score_csv(committed_file)
    legacy_table <- read_legacy_z_score_csv(legacy_file)

    committed_difference <- compare_z_score_tables(
      regenerated_table,
      committed_table,
      paste(target$cohort_label, target$pathway_id, "regenerated-versus-committed"),
      tolerance
    )
    legacy_difference <- compare_z_score_tables(
      regenerated_table,
      legacy_table,
      paste(target$cohort_label, target$pathway_id, "regenerated-versus-legacy"),
      tolerance
    )

    z_values <- as.matrix(regenerated_table[, ko_z_score_treatment_levels, drop = FALSE])
    storage.mode(z_values) <- "double"
    complete_rows <- stats::complete.cases(z_values)
    all_missing_rows <- rowSums(is.na(z_values)) == ncol(z_values)
    check_ko_z_score(
      all(complete_rows | all_missing_rows),
      paste(target$pathway_id, "contains a partially missing z-score row.")
    )
    check_ko_z_score(
      all(abs(rowMeans(z_values[complete_rows, , drop = FALSE])) <= tolerance),
      paste(target$pathway_id, "has a z-score row whose mean is not zero.")
    )
    check_ko_z_score(
      all(abs(apply(z_values[complete_rows, , drop = FALSE], 1L, stats::sd) - 1) <= tolerance),
      paste(target$pathway_id, "has a z-score row whose sample SD is not one.")
    )

    data.frame(
      cohort = target$cohort_label,
      pathway_id = target$pathway_id,
      ko_rows = nrow(regenerated_table),
      undefined_constant_rows = sum(all_missing_rows),
      regenerated_vs_committed_max_abs_diff = committed_difference,
      regenerated_vs_legacy_max_abs_diff = legacy_difference,
      stringsAsFactors = FALSE
    )
  })

  summary <- do.call(rbind, comparisons)
  rownames(summary) <- NULL
  check_ko_z_score(nrow(summary) == 19L, "The validation did not cover all 19 source tables.")
  message(
    "Validated ", nrow(summary), " KEGG-pathway KO z-score tables (",
    sum(summary$ko_rows), " pathway-KO rows); maximum legacy difference = ",
    format(max(summary$regenerated_vs_legacy_max_abs_diff), scientific = TRUE),
    "; canonical-versus-KO2 input difference = ",
    format(max(canonical_input_comparison$max_absolute_difference), scientific = TRUE), "."
  )
  attr(summary, "canonical_input_comparison") <- canonical_input_comparison
  summary
}

if (sys.nframe() == 0L) {
  validation <- validate_kegg_pathway_ko_z_scores()
  print(validation, row.names = FALSE)
}
