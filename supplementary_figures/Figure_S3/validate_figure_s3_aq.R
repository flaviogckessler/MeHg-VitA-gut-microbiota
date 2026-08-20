# Validation for the Figure S3 panels A-Q workflow.

# Load Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
})

find_figure_s3_generator <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "generate_figure_s3_aq.R"),
      file.path(
        path,
        "supplementary_figures",
        "Figure_S3",
        "generate_figure_s3_aq.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_s3_aq.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s3_generator())
validate_supplementary_abundance_figure(
  figure_s3_specification,
  generate_figure_s3_aq
)
message(
  "Figure S3 panels A-Q passed source, cohort, nonparametric-statistics, ",
  "annotation, and PNG-rendering checks."
)
