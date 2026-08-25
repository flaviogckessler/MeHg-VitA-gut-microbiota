# Reproducible KO-only generator for Figure 6 panel B (F1, PND365).
# The historical reference is scripts/legacy_scripts/picrust_legacy.R;
# it is retained for provenance and is not sourced at runtime.

# Load Packages ----
suppressPackageStartupMessages({
  library(ggpicrust2)
  library(ALDEx2)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

find_figure_6_directory <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(path, file.path(path, "figures", "Figure_6"))
    for (candidate in candidates) {
      common <- file.path(candidate, "figure_6_common.R")
      if (file.exists(common)) return(candidate)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate figures/Figure_6/figure_6_common.R.", call. = FALSE)
    }
    path <- parent
  }
}

figure_6_directory <- find_figure_6_directory()
source(file.path(figure_6_directory, "figure_6_common.R"))

prepare_figure_6_b_data <- function(
    project_root = find_project_root(),
    seed = 20260809L) {
  analysis <- prepare_figure_6_analysis(
    project_root = project_root,
    batch_id = "batch02_illumina-miseq_f1-f2_pnd365",
    generation = "F1",
    seed = seed
  )
  threshold <- 0.10
  # ko00965 (Betalain biosynthesis) crosses the relaxed boundary in the
  # seeded rerun but was not retained in the legacy-unseeded final panel.
  excluded_legacy_pathways <- c("ko00514", "ko00944", "ko00965")
  expected_features <- "ko03050" # Proteasome
  threshold_selected <- analysis$daa[
    analysis$daa$p_adjust <= threshold &
      !analysis$daa$feature %in% excluded_legacy_pathways,
    ,
    drop = FALSE
  ]
  unexpected <- setdiff(threshold_selected$feature, expected_features)
  selected <- analysis$daa[
    match(expected_features, analysis$daa$feature),
    ,
    drop = FALSE
  ]
  if (anyNA(selected$feature) || length(unexpected) > 0L ||
      any(selected$p_adjust > threshold + 0.01)) {
    stop(
      paste0(
        "Panel B pathway selection differs from the final-figure specification. ",
        "Threshold-selected: ", paste(threshold_selected$feature, collapse = ", "),
        "; expected author-defined sixfold omnibus values: ",
        paste(format_figure_6_p_value(selected$p_adjust), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  analysis$selected <- selected
  analysis$expected_features <- expected_features
  analysis$threshold <- threshold
  analysis$threshold_tolerance <- 0.01
  analysis$excluded_legacy_pathways <- excluded_legacy_pathways
  analysis
}

generate_figure_6_b <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_6", "output"),
    seed = 20260809L) {
  analysis <- prepare_figure_6_b_data(seed = seed)
  group <- analysis$input$metadata$treatment[
    match(colnames(analysis$abundance), analysis$input$metadata$sample_id)
  ]
  plot <- pathway_errorbar3(
    abundance = analysis$abundance,
    daa_results_df = analysis$daa,
    Group = group,
    ko_to_kegg = TRUE,
    p_values_threshold = analysis$threshold,
    order = "pathway_class",
    select = analysis$expected_features,
    p_value_bar = FALSE,
    colors = figure_6_treatment_colors,
    x_lab = "pathway_name",
    panel_letter = "B",
    title = "F1 -PND365",
    show_legend = FALSE
  )
  file <- file.path(output_dir, "Figure6B_generated.png")
  save_figure_6_panel(plot, file, width_px = 2048, height_px = 265)
  invisible(list(
    file = file,
    plot = plot,
    analysis = analysis,
    selected_daa = attr(plot, "selected_daa"),
    relative_abundance_summary = attr(plot, "relative_abundance_summary")
  ))
}

if (sys.nframe() == 0L) {
  result <- generate_figure_6_b()
  message("Generated Figure 6 panel B: ", normalizePath(result$file, winslash = "/"))
  print(
    result$selected_daa[, c(
      "feature", "pathway_name", "raw_omnibus_p", "package_p_adjust",
      "author_sixfold_omnibus_p"
    )],
    row.names = FALSE
  )
}
