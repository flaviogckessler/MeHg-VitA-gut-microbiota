# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(rlang)
  library(tools)
})

# - Load the figure workflow
find_figure_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidates <- c(
      file.path(path, "generate_figure_3_hk.R"),
      file.path(path, "figures", "Figure_3", "generate_figure_3_hk.R")
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) return(existing[[1]])

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_3_hk.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_script())

# - Validation helpers
check_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

# - Input validation
project_root <- find_project_root()
abundance_data <- prepare_abundance_data(project_root)
panel_data <- prepare_panel_hk_data(abundance_data)
expected_counts <- c(H = 18L, I = 34L, J = 24L, K = 24L)
observed_counts <- vapply(
  panel_data,
  function(data) n_distinct(data$sample_id),
  integer(1)
)

check_true(
  identical(observed_counts, expected_counts),
  "Panels H-K do not contain the expected sample counts."
)
check_true(
  all(vapply(panel_data, function(data) {
    identical(levels(data$phylum), phylum_levels)
  }, logical(1))),
  "The four phyla are not in the required order."
)
check_true(
  all(vapply(panel_data, function(data) {
    nrow(data) == n_distinct(data$sample_id) * length(phylum_levels)
  }, logical(1))),
  "A panel does not contain one value per sample and phylum."
)
check_true(
  all(vapply(panel_data, function(data) {
    all(is.finite(data$relative_abundance)) &&
      all(data$relative_abundance >= 0) &&
      all(data$relative_abundance <= 1)
  }, logical(1))),
  "One or more relative-abundance values are invalid."
)
check_true(
  !any(c("10F", "22F") %in% abundance_data$batch1$sample_id),
  "Excluded batch 1 samples must not be analyzed."
)

# - Statistical validation
statistics <- map(panel_data, calculate_panel_hk_statistics)
expected_methods <- list(
  H = c("parametric", "parametric", "nonparametric", "nonparametric"),
  I = rep("nonparametric", 4),
  J = c("nonparametric", "parametric", "nonparametric", "nonparametric"),
  K = rep("nonparametric", 4)
)

check_true(
  all(map2_lgl(statistics, expected_methods, function(result, expected) {
    identical(result$shapiro$method, expected)
  })),
  "Shapiro-Wilk classifications do not match the expected test methods."
)
check_true(
  all(vapply(statistics, function(result) {
    nrow(result$treatment) == 24 && all(is.finite(result$treatment$p.adj))
  }, logical(1))),
  "One or more treatment multiple-comparison tests failed."
)
check_true(
  is.null(statistics$H$sex),
  "Panel H must not run a between-sex comparison."
)
check_true(
  all(vapply(statistics[c("I", "J", "K")], function(result) {
    nrow(result$sex) == 4 && all(is.finite(result$sex$p))
  }, logical(1))),
  "One or more between-sex tests failed."
)

# - Output validation
validation_dir <- file.path(tempdir(), "figure_3_hk_validation")
result <- generate_figure_3_hk(validation_dir)
expected_files <- unname(result$files)
expected_names <- paste0("Figure3", LETTERS[8:11], "_generated.png")

check_true(all(file.exists(expected_files)), "One or more PNG panels were not generated.")
check_true(
  identical(basename(expected_files), expected_names),
  "Generated panel filenames do not match the expected names."
)
check_true(
  all(tolower(tools::file_ext(expected_files)) == "png"),
  "Generated figure files must use PNG format."
)
check_true(
  all(file.info(expected_files)$size > 1000 & file.info(expected_files)$size < 500000),
  "One or more generated PNG panels are empty or unnecessarily large."
)
check_true(
  all(vapply(result$plots, inherits, logical(1), what = "ggplot")),
  "One or more panel objects are not ggplot objects."
)
check_true(
  all(vapply(result$plots, function(plot) {
    identical(rlang::as_label(plot$mapping$x), "treatment")
  }, logical(1))),
  "Each phylum facet must distribute boxplots across treatment positions."
)

message("Figure 3 panels H-K passed all validation checks.")
