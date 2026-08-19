# Validation for the Figure S4 panels A-AN workflow.

find_figure_s4_generator <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "generate_figure_s4_a_an.R"),
      file.path(
        path,
        "supplementary_figures",
        "Figure_S4",
        "generate_figure_s4_a_an.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_s4_a_an.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s4_generator())
validate_supplementary_abundance_figure(
  figure_s4_specification,
  generate_figure_s4_a_an
)
message(
  "Figure S4 panels A-AN passed source, cohort, nonparametric-statistics, ",
  "annotation, and PNG-rendering checks."
)
