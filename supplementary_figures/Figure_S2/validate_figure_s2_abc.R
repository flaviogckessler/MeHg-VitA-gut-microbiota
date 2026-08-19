# Validation for the Figure S2 panels A-C workflow.

find_figure_s2_generator <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "generate_figure_s2_abc.R"),
      file.path(
        path,
        "supplementary_figures",
        "Figure_S2",
        "generate_figure_s2_abc.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_s2_abc.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s2_generator())
validate_supplementary_abundance_figure(
  figure_s2_specification,
  generate_figure_s2_abc
)
message(
  "Figure S2 panels A-C passed source, cohort, nonparametric-statistics, ",
  "annotation, and PNG-rendering checks."
)
