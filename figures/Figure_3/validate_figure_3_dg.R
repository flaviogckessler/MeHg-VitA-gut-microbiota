# - Load the figure workflow
find_figure_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidates <- c(
      file.path(path, "generate_figure_3_dg.R"),
      file.path(path, "figures", "Figure_3", "generate_figure_3_dg.R")
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0) return(existing[[1]])

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_3_dg.R.", call. = FALSE)
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
alpha_data <- prepare_alpha_data(project_root)
panel_data <- prepare_panel_dg_data(alpha_data)
expected_counts <- c(D = 18L, E = 34L, F = 24L, G = 24L)
observed_counts <- vapply(
  panel_data,
  function(data) n_distinct(data$sample_id),
  integer(1)
)

check_true(
  identical(observed_counts, expected_counts),
  "Panels D-G do not contain the expected sample counts."
)
check_true(
  !"t3" %in% alpha_data$sample_id,
  "The undocumented batch 2 sample t3 must not be analyzed."
)
check_true(
  all(vapply(panel_data, function(data) all(is.finite(data$Shannon)), logical(1))),
  "One or more panels contain invalid Shannon values."
)
check_true(
  identical(as.character(unique(panel_data$D$sex)), "Female"),
  "Panel D must contain F0 dams only."
)
check_true(
  all(vapply(panel_data[c("E", "F", "G")], function(data) {
    n_distinct(data$sex) == 2
  }, logical(1))),
  "Panels E-G must contain both sexes."
)

# - Statistical validation
methods <- c(
  D = "parametric",
  E = "nonparametric",
  F = "nonparametric",
  G = "parametric"
)
statistics <- map2(panel_data, methods, calculate_panel_statistics)
check_true(
  is.null(statistics$D$sex),
  "Panel D must not run a between-sex comparison."
)
check_true(
  all(vapply(statistics[c("E", "F", "G")], function(result) {
    nrow(result$sex) == 1 && is.finite(result$sex$p)
  }, logical(1))),
  "One or more between-sex tests failed."
)

# - Output validation
validation_dir <- file.path(tempdir(), "figure_3_dg_validation")
result <- generate_figure_3_dg(validation_dir)
expected_files <- unname(result$files)
expected_names <- paste0("Figure3", LETTERS[4:7], "_generated.png")

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

message("Figure 3 panels D-G passed all validation checks.")
