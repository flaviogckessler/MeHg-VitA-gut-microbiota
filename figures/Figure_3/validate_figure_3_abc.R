# Load Packages ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(tools)
})

# - Load the figure workflow
find_figure_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidates <- c(
      file.path(path, "generate_figure_3_abc.R"),
      file.path(path, "figures", "Figure_3", "generate_figure_3_abc.R")
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) return(existing[[1]])

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_3_abc.R.", call. = FALSE)
    }
    path <- parent
  }
}

script_path <- find_figure_script()
source(script_path)

# - Validation helpers
check_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

# - Input validation
project_root <- find_project_root()
panel_a_data <- prepare_panel_a_data(project_root)
panels_bc_data <- prepare_panels_bc_data(project_root)

check_true(
  n_distinct(panel_a_data$sample_id) == 52,
  "Panel A must contain 52 included samples."
)
check_true(
  !any(c("10F", "22F") %in% panel_a_data$sample_id),
  "Panel A contains samples excluded for low sequencing quality."
)
check_true(
  all(count(panel_a_data, sample_id)$n > 0),
  "Panel A contains a sample without abundance records."
)
check_true(
  all(count(panels_bc_data, generation)$n > 0),
  "Panels B or C lack abundance records."
)
check_true(
  identical(
    count(distinct(panels_bc_data, sample_id, generation), generation)$n,
    c(24L, 24L)
  ),
  "Panels B and C must each contain 24 samples."
)

sample_totals <- panels_bc_data %>%
  group_by(sample_id) %>%
  summarise(total = sum(abundance), .groups = "drop")
check_true(
  all(abs(sample_totals$total - 1) < 1e-8),
  "Batch 2 relative abundances must sum to one for every sample."
)

# - Output validation
validation_dir <- file.path(tempdir(), "figure_3_abc_validation")
result <- generate_figure_3_abc(validation_dir)
expected_files <- unname(result$files)
expected_names <- c(
  "Figure3A_generated.png",
  "Figure3B_generated.png",
  "Figure3C_generated.png"
)

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
  all(file.info(expected_files)$size > 1000),
  "One or more generated PNG panels are empty."
)
check_true(
  all(vapply(result$plots, inherits, logical(1), what = "ggplot")),
  "One or more panel objects are not ggplot objects."
)

message("Figure 3 panels A-C passed all validation checks.")
