# End-to-end validation for the portable Figure S8-S16 heatmap renderers.

# Load Packages ----
suppressPackageStartupMessages({
  library(graphics)
  library(grDevices)
  library(stats)
  library(utils)
})

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

validate_figure_s8_s16 <- function(
    project_root = find_supplementary_heatmap_project_root()) {
  temporary_root <- tempfile("figure_s8_s16_validation_")
  dir.create(temporary_root, recursive = TRUE, showWarnings = FALSE)

  summaries <- list()
  summary_index <- 0L
  for (figure_id in names(supplementary_heatmap_configurations)) {
    figure_configuration <- supplementary_heatmap_configurations[[figure_id]]
    wrapper_pattern <- if (identical(figure_id, "S15")) {
      "generate_figure_s15_ij.R"
    } else {
      paste0("generate_figure_", tolower(figure_id), "_e.R")
    }
    wrapper_file <- file.path(
      project_root, "supplementary_figures", paste0("Figure_", figure_id),
      wrapper_pattern
    )
    if (!file.exists(wrapper_file)) {
      stop("Missing per-figure generator: ", wrapper_file, call. = FALSE)
    }
    tryCatch(
      parse(file = wrapper_file),
      error = function(error) stop("Could not parse ", wrapper_file, ": ", error$message, call. = FALSE)
    )

    temporary_output <- file.path(temporary_root, figure_id)
    regenerated <- generate_supplementary_heatmap_figure(
      figure_id = figure_id,
      output_dir = temporary_output,
      project_root = project_root
    )
    for (panel_name in names(regenerated)) {
      result <- regenerated[[panel_name]]
      panel_configuration <- figure_configuration$panels[[panel_name]]
      committed_file <- supplementary_heatmap_output_file(
        project_root, figure_id, panel_configuration
      )
      if (!file.exists(committed_file)) {
        stop("Missing committed generated panel: ", committed_file, call. = FALSE)
      }
      if (!identical(
        unname(read_supplementary_png_dimensions(committed_file)),
        c(panel_configuration$width, panel_configuration$height)
      )) {
        stop("Committed panel dimensions changed: ", committed_file, call. = FALSE)
      }

      summary_index <- summary_index + 1L
      summaries[[summary_index]] <- data.frame(
        figure = figure_id,
        panel = panel_name,
        cohort = panel_configuration$cohort_label,
        pathway_id = panel_configuration$pathway_id,
        source_ko_rows = panel_configuration$expected_source_ko_count,
        rendered_ko_rows = length(panel_configuration$reference_row_order),
        omitted_nonfinite_rows = length(panel_configuration$omitted_nonfinite_ko_ids),
        width = panel_configuration$width,
        height = panel_configuration$height,
        stringsAsFactors = FALSE
      )
    }
  }

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  if (nrow(summary) != 10L) {
    stop("Expected validation coverage for 10 heatmap panels.", call. = FALSE)
  }
  message(
    "Validated ", nrow(summary), " heatmap panels across Figures S8-S16 (",
    sum(summary$rendered_ko_rows), " finite KO rows; ",
    sum(summary$omitted_nonfinite_rows), " explicitly omitted undefined rows)."
  )
  summary
}

if (sys.nframe() == 0L) {
  print(validate_figure_s8_s16(), row.names = FALSE)
}
