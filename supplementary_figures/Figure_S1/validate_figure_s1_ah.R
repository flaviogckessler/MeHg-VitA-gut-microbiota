# Validation for the portable Figure S1 panels A-H workflow.

# Load Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(tools)
})

find_figure_s1_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidates <- c(
      file.path(path, "generate_figure_s1_ah.R"),
      file.path(
        path,
        "supplementary_figures",
        "Figure_S1",
        "generate_figure_s1_ah.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1]])

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_s1_ah.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s1_script())

check_figure_s1 <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_figure_s1_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)

  signature <- readBin(connection, what = "raw", n = 8L)
  expected_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (!identical(signature, expected_signature)) {
    stop("Figure S1 output is not a valid PNG: ", filename, call. = FALSE)
  }

  chunk_length <- readBin(
    connection,
    what = "integer",
    n = 1L,
    size = 4L,
    endian = "big"
  )
  chunk_type <- rawToChar(readBin(connection, what = "raw", n = 4L))
  if (!identical(chunk_length, 13L) || !identical(chunk_type, "IHDR")) {
    stop("Figure S1 PNG does not begin with an IHDR chunk.", call. = FALSE)
  }

  c(
    width = readBin(
      connection,
      what = "integer",
      n = 1L,
      size = 4L,
      endian = "big"
    ),
    height = readBin(
      connection,
      what = "integer",
      n = 1L,
      size = 4L,
      endian = "big"
    )
  )
}

project_root <- find_figure_s1_project_root()
source_files <- figure_s1_source_files(project_root)
check_figure_s1(
  all(file.exists(unlist(source_files, use.names = FALSE))),
  "One or more canonical Figure S1 sources are missing."
)

alpha_data <- prepare_figure_s1_alpha_data(project_root)
check_figure_s1(
  nrow(alpha_data) == 100L && !anyDuplicated(alpha_data$sample_id),
  "Figure S1 must use 100 unique primary-analysis samples."
)
check_figure_s1(
  identical(
    as.integer(table(factor(
      alpha_data$batch_id,
      levels = unname(figure_s1_batch_ids)
    ))),
    c(52L, 48L)
  ),
  "Figure S1 batch sample counts are not 52 and 48."
)
check_figure_s1(
  !any(c("10F", "22F", "t3") %in% alpha_data$sample_id),
  "Excluded or undocumented samples entered the Figure S1 analysis."
)

unmatched_ids <- attr(alpha_data, "source_ids_without_manifest")
excluded_ids <- attr(alpha_data, "excluded_nonprimary_ids")
check_figure_s1(
  identical(unmatched_ids[[figure_s1_batch_ids[["batch1"]]]], character(0)) &&
    identical(unmatched_ids[[figure_s1_batch_ids[["batch2"]]]], "t3"),
  "The expected undocumented batch 02 alpha-diversity ID is not exactly t3."
)
check_figure_s1(
  identical(
    excluded_ids[[figure_s1_batch_ids[["batch1"]]]],
    c("10F", "22F")
  ) &&
    identical(
      excluded_ids[[figure_s1_batch_ids[["batch2"]]]],
      character(0)
    ),
  "Primary-analysis exclusions no longer match 10F and 22F."
)
check_figure_s1(
  all(is.finite(alpha_data$Observed)) &&
    all(is.finite(alpha_data$InvSimpson)) &&
    all(alpha_data$Observed > 0) &&
    all(alpha_data$InvSimpson >= 1),
  "Figure S1 contains invalid Observed or InvSimpson values."
)

panel_data <- prepare_figure_s1_panel_data(alpha_data)
expected_counts <- c(A = 18L, B = 34L, C = 24L, D = 24L,
                     E = 18L, F = 34L, G = 24L, H = 24L)
observed_counts <- vapply(
  panel_data,
  function(data) length(unique(data$sample_id)),
  integer(1)
)
check_figure_s1(
  identical(observed_counts, expected_counts),
  "Figure S1 panels do not contain the expected sample counts."
)
check_figure_s1(
  all(vapply(panel_data[c("A", "E")], function(data) {
    identical(as.character(unique(data$sex)), "Female")
  }, logical(1))),
  "Figure S1 dams panels must contain females only."
)
check_figure_s1(
  all(vapply(panel_data[setdiff(names(panel_data), c("A", "E"))], function(data) {
    identical(sort(as.character(unique(data$sex))), figure_s1_sex_levels)
  }, logical(1))),
  "Figure S1 offspring panels must contain both sexes."
)
check_figure_s1(
  all(vapply(panel_data, function(data) {
    identical(
      sort(as.character(unique(data$treatment))),
      sort(figure_s1_treatment_levels)
    )
  }, logical(1))),
  "Every Figure S1 panel must contain all four treatment groups."
)

expected_metrics <- c(
  A = "Observed", B = "Observed", C = "Observed", D = "Observed",
  E = "InvSimpson", F = "InvSimpson", G = "InvSimpson", H = "InvSimpson"
)
expected_methods <- c(
  A = "parametric", B = "nonparametric", C = "nonparametric", D = "parametric",
  E = "parametric", F = "parametric", G = "parametric", H = "parametric"
)
observed_metrics <- vapply(
  figure_s1_panel_configurations,
  `[[`,
  character(1),
  "metric"
)
observed_methods <- vapply(
  figure_s1_panel_configurations,
  `[[`,
  character(1),
  "method"
)
check_figure_s1(
  identical(observed_metrics, expected_metrics),
  "Figure S1 panel-to-metric mapping changed."
)
check_figure_s1(
  identical(observed_methods, expected_methods),
  "Figure S1 panel statistical methods changed."
)

statistics <- Map(function(data, configuration) {
  calculate_figure_s1_statistics(
    data,
    configuration$metric,
    configuration$method
  )
}, panel_data, figure_s1_panel_configurations)
check_figure_s1(
  all(vapply(statistics[c("A", "E")], function(result) {
    is.null(result$sex)
  }, logical(1))),
  "Dams panels must not run a between-sex test."
)
check_figure_s1(
  all(vapply(statistics[setdiff(names(statistics), c("A", "E"))], function(result) {
    nrow(result$sex) == 1L && is.finite(result$sex$p) &&
      result$sex$p >= 0 && result$sex$p <= 1
  }, logical(1))),
  "One or more Figure S1 between-sex tests failed."
)
expected_sex_statistics <- c(
  B = 169,
  C = 68,
  D = 1.4476089,
  F = -0.5100575,
  G = 0.9252663,
  H = 2.5634399
)
expected_sex_p_values <- c(
  B = 0.3978182761390764,
  C = 0.8398256197101182,
  D = 0.1623679533760797,
  F = 0.6135141297378561,
  G = 0.3651602176666696,
  H = 0.01888787291363542
)
observed_sex_statistics <- vapply(
  statistics[names(expected_sex_statistics)],
  function(result) result$sex$statistic,
  numeric(1)
)
observed_sex_p_values <- vapply(
  statistics[names(expected_sex_p_values)],
  function(result) result$sex$p,
  numeric(1)
)
check_figure_s1(
  max(abs(observed_sex_statistics - expected_sex_statistics)) < 1e-6 &&
    max(abs(observed_sex_p_values - expected_sex_p_values)) < 1e-12,
  "Figure S1 between-sex statistical checkpoints changed."
)
expected_treatment_test_counts <- c(
  A = 6L, B = 12L, C = 12L, D = 12L,
  E = 6L, F = 12L, G = 12L, H = 12L
)
observed_treatment_test_counts <- vapply(
  statistics,
  function(result) nrow(result$treatment),
  integer(1)
)
check_figure_s1(
  identical(observed_treatment_test_counts, expected_treatment_test_counts),
  "Figure S1 did not evaluate all expected within-sex treatment pairs."
)
check_figure_s1(
  all(vapply(statistics, function(result) {
    all(is.finite(result$treatment$p.adj)) &&
      all(result$treatment$p.adj >= 0) &&
      all(result$treatment$p.adj <= 1)
  }, logical(1))),
  "One or more Figure S1 adjusted treatment p-values are invalid."
)
expected_significant_treatment_keys <- list(
  A = character(0),
  B = character(0),
  C = character(0),
  D = character(0),
  E = character(0),
  F = c(
    "Female|CN|VitA",
    "Female|CN|MeHg",
    "Female|CN|VitA+MeHg"
  ),
  G = "Female|VitA|MeHg",
  H = character(0)
)
observed_significant_treatment_keys <- lapply(statistics, function(result) {
  significant <- result$treatment[
    result$treatment$p.adj < 0.05,
    ,
    drop = FALSE
  ]
  paste(significant$sex, significant$group1, significant$group2, sep = "|")
})
check_figure_s1(
  identical(
    observed_significant_treatment_keys,
    expected_significant_treatment_keys
  ),
  "Figure S1 significant within-sex treatment comparisons changed."
)

validation_dir <- file.path(tempdir(), "figure_s1_ah_validation")
result <- generate_figure_s1_ah(validation_dir, project_root)
expected_names <- paste0("FigureS1", LETTERS[1:8], "_generated.png")
check_figure_s1(
  identical(basename(unname(result$files)), expected_names),
  "Generated Figure S1 panel filenames do not match the repository convention."
)
check_figure_s1(
  all(file.exists(result$files)) && all(file.info(result$files)$size > 1000),
  "One or more Figure S1 validation PNGs were not generated correctly."
)
check_figure_s1(
  all(vapply(result$plots, inherits, logical(1), what = "ggplot")),
  "One or more Figure S1 panel objects are not ggplot objects."
)

expected_dimensions <- lapply(
  c(600L, 1200L, 1200L, 1200L, 600L, 1200L, 1200L, 1200L),
  function(width) c(width = width, height = 750L)
)
validation_dimensions <- lapply(
  unname(result$files),
  read_figure_s1_png_dimensions
)
check_figure_s1(
  identical(validation_dimensions, expected_dimensions),
  "One or more generated Figure S1 panels have unexpected dimensions."
)

official_files <- file.path(
  project_root,
  "supplementary_figures",
  "Figure_S1",
  "output",
  expected_names
)
check_figure_s1(
  all(file.exists(official_files)) && all(file.info(official_files)$size > 1000),
  "One or more committed Figure S1 panel PNGs are missing or empty."
)
official_dimensions <- lapply(official_files, read_figure_s1_png_dimensions)
check_figure_s1(
  identical(official_dimensions, expected_dimensions),
  "One or more committed Figure S1 panels have unexpected dimensions."
)

reference_composite <- file.path(
  project_root,
  "supplementary_figures",
  "Figure_S1",
  "output",
  "Figure_S1_FinalFigure.png"
)
check_figure_s1(
  file.exists(reference_composite) && file.info(reference_composite)$size > 1000,
  "The supplied Figure S1 reference composite is missing or empty."
)
check_figure_s1(
  identical(
    read_figure_s1_png_dimensions(reference_composite),
    c(width = 2650L, height = 1145L)
  ),
  "The supplied Figure S1 reference-composite dimensions changed."
)
check_figure_s1(
  !any(tolower(tools::file_ext(c(result$files, official_files))) %in% c("tif", "tiff")),
  "Figure S1 generation must not create TIFF files."
)

message("Figure S1 panels A-H passed all validation checks.")
