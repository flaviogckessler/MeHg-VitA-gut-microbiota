# Reproducible Figure S14 panel E: F1-PND21/map00930 KO z-score heatmap.
find_common_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidate <- file.path(path, "supplementary_figures", "figure_s8_s16_common.R")
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not locate figure_s8_s16_common.R.", call. = FALSE)
    path <- parent
  }
}
source(find_common_script())
generate_figure_s14_e <- function(output_dir = NULL) {
  generate_supplementary_heatmap_figure("S14", output_dir)
}
if (sys.nframe() == 0L) {
  result <- generate_figure_s14_e()
  message("Generated Figure S14 panel E: ", normalizePath(result$E$file, winslash = "/"))
}
