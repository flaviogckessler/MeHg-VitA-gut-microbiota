# Validation for the Figure S5 panels A-O workflow.

find_figure_s5_generator <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "generate_figure_s5_ao.R"),
      file.path(
        path,
        "supplementary_figures",
        "Figure_S5",
        "generate_figure_s5_ao.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_s5_ao.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s5_generator())
validate_supplementary_abundance_figure(
  figure_s5_specification,
  generate_figure_s5_ao
)
message(
  "Figure S5 panels A-O passed source, cohort, nonparametric-statistics, ",
  "annotation, and PNG-rendering checks."
)
