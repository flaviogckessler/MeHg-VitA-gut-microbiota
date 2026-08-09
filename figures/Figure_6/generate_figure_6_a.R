# Reproducible KO-only generator for Figure 6 panel A (F1, PND21).
# The historical reference is scripts/legacy_scripts/Picrust_mestrado_Alberto.R;
# it is retained for provenance and is not sourced at runtime.

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

prepare_figure_6_a_data <- function(
    project_root = find_project_root(),
    seed = 20260809L) {
  analysis <- prepare_figure_6_analysis(
    project_root = project_root,
    batch_id = "batch01_iontorrent_dams-f0_f1-pnd21",
    generation = "F1",
    seed = seed
  )
  threshold <- 0.05
  excluded_legacy_pathways <- c(
    "ko00514", "ko00944", "ko04011", "ko04962", "ko05219"
  )
  expected_features <- c(
    "ko03050", # Proteasome
    "ko03450", # Non-homologous end-joining
    "ko00310", # Lysine degradation
    "ko00190", # Oxidative phosphorylation
    "ko00740", # Riboflavin metabolism
    "ko00760", # Nicotinate and nicotinamide metabolism
    "ko00430", # Taurine and hypotaurine metabolism
    "ko00480", # Glutathione metabolism
    "ko00623", # Toluene degradation
    "ko00930"  # Caprolactam degradation
  )
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
        "Panel A pathway selection differs from the final-figure specification. ",
        "Threshold-selected: ", paste(threshold_selected$feature, collapse = ", "),
        "; expected adjusted p-values: ",
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

generate_figure_6_a <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_6", "output"),
    seed = 20260809L) {
  analysis <- prepare_figure_6_a_data(seed = seed)
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
    panel_letter = "A",
    title = "F1 -PND21",
    show_legend = TRUE
  )
  file <- file.path(output_dir, "Figure6A_generated.png")
  save_figure_6_panel(plot, file, width_px = 2048, height_px = 760)
  invisible(list(
    file = file,
    plot = plot,
    analysis = analysis,
    selected_daa = attr(plot, "selected_daa"),
    relative_abundance_summary = attr(plot, "relative_abundance_summary")
  ))
}

if (sys.nframe() == 0L) {
  result <- generate_figure_6_a()
  message("Generated Figure 6 panel A: ", normalizePath(result$file, winslash = "/"))
  print(
    result$selected_daa[, c("feature", "pathway_name", "p_values", "p_adjust")],
    row.names = FALSE
  )
}
