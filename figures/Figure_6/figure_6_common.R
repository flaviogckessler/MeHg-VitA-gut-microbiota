# Shared KO-only analysis and plotting helpers for Figure 6 panels A-C.
#
# This file is sourced by the three panel generators. It replaces the fragile
# multi-group bookkeeping in the legacy pathway_errorbar3() while preserving
# the analysis choices used for the final figure:
#   * PICRUSt2 KO predictions aggregated to KEGG pathways with ggpicrust2;
#   * one four-treatment omnibus ALDEx2 GLM test per pathway;
#   * the author's additional conservative sixfold rule applied to each raw
#     omnibus p-value, with no post-hoc pairwise testing; and
#   * per-sample pathway relative abundance summarized as mean +/- SD.

# Load Packages ----
suppressPackageStartupMessages({
  library(ggpicrust2)
  library(ALDEx2)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stats)
  library(utils)
})

figure_6_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
figure_6_treatment_colors <- c(
  CN = "#43DFE8",
  VitA = "#80E500",
  MeHg = "#F51F26",
  `VitA+MeHg` = "#454545"
)

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata", "sample_manifest_all_batches.csv"),
      file.path(path, "results", "functional_prediction"),
      file.path(path, "figures", "Figure_6")
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the New_repo project root.", call. = FALSE)
    }
    path <- parent
  }
}

require_figure_6_packages <- function() {
  packages <- c("ggpicrust2", "ALDEx2", "ggplot2", "patchwork", "scales")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      paste(
        "Figure 6 requires the following installed R packages:",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(packages)
}

is_primary_analysis_sample <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

read_figure_6_input <- function(project_root, batch_id, generation) {
  abundance_file <- file.path(
    project_root,
    "results", "functional_prediction", batch_id,
    "kegg_ortholog_abundance.tsv"
  )
  metadata_file <- file.path(
    project_root, "data", "metadata", "sample_manifest_all_batches.csv"
  )

  if (!file.exists(abundance_file)) {
    stop("Missing canonical KO abundance file: ", abundance_file, call. = FALSE)
  }
  if (!file.exists(metadata_file)) {
    stop("Missing canonical sample manifest: ", metadata_file, call. = FALSE)
  }

  abundance <- utils::read.delim(
    abundance_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  metadata <- utils::read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )

  required_metadata <- c(
    "sample_id", "batch_id", "generation", "treatment",
    "include_primary_analysis"
  )
  if (!all(required_metadata %in% names(metadata))) {
    stop("The canonical manifest is missing required Figure 6 columns.", call. = FALSE)
  }
  if (ncol(abundance) < 2L || names(abundance)[1L] != "ID") {
    stop("The KO abundance table must start with an ID column.", call. = FALSE)
  }

  metadata <- metadata[
    metadata$batch_id == batch_id &
      metadata$generation == generation &
      is_primary_analysis_sample(metadata$include_primary_analysis),
    ,
    drop = FALSE
  ]
  metadata$treatment <- factor(
    metadata$treatment,
    levels = figure_6_treatment_levels
  )
  if (nrow(metadata) == 0L || anyNA(metadata$treatment)) {
    stop("No valid samples/treatment labels remain for this panel.", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("The canonical manifest contains duplicated sample IDs.", call. = FALSE)
  }

  missing_samples <- setdiff(metadata$sample_id, names(abundance))
  if (length(missing_samples) > 0L) {
    stop(
      "Included samples missing from KO abundance table: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  sample_ids <- metadata$sample_id
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
  ko <- abundance[, c("ID", sample_ids), drop = FALSE]
  ko$ID <- sub("^ko:", "", ko$ID)
  if (any(!grepl("^K[0-9]{5}$", ko$ID))) {
    stop("The normalized KO table contains malformed KEGG Ortholog IDs.", call. = FALSE)
  }
  if (anyDuplicated(ko$ID)) {
    stop("The normalized KO table contains duplicated KO IDs.", call. = FALSE)
  }
  numeric_matrix <- as.matrix(ko[, sample_ids, drop = FALSE])
  storage.mode(numeric_matrix) <- "double"
  if (any(!is.finite(numeric_matrix)) || any(numeric_matrix < 0)) {
    stop("The KO abundance table contains invalid numeric values.", call. = FALSE)
  }
  ko[, sample_ids] <- numeric_matrix

  list(
    ko = ko,
    metadata = metadata,
    sample_ids = sample_ids,
    source_files = c(abundance = abundance_file, metadata = metadata_file)
  )
}

aggregate_ko_to_kegg <- function(input) {
  require_figure_6_packages()
  invisible(utils::capture.output({
    abundance <- suppressWarnings(
      suppressMessages(ggpicrust2::ko2kegg_abundance(data = input$ko))
    )
  }))
  abundance <- as.data.frame(abundance, check.names = FALSE)
  abundance <- abundance[, input$sample_ids, drop = FALSE]
  values <- as.matrix(abundance)
  storage.mode(values) <- "double"
  if (nrow(values) == 0L || any(!is.finite(values))) {
    stop("KO-to-KEGG aggregation produced invalid pathway abundances.", call. = FALSE)
  }
  abundance[] <- values
  abundance
}

run_figure_6_daa <- function(kegg_abundance, metadata, seed = 20260809L) {
  require_figure_6_packages()
  metadata <- metadata[match(colnames(kegg_abundance), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id)) {
    stop("KEGG abundance columns could not be aligned to metadata.", call. = FALSE)
  }
  metadata$treatment <- factor(
    metadata$treatment,
    levels = figure_6_treatment_levels
  )

  set.seed(seed)
  daa <- suppressMessages(
    ggpicrust2::pathway_daa(
      abundance = kegg_abundance,
      metadata = metadata[, c("sample_id", "treatment"), drop = FALSE],
      group = "treatment",
      daa_method = "ALDEx2",
      select = NULL,
      reference = NULL
    )
  )
  daa <- daa[daa$method == "ALDEx2_glm test", , drop = FALSE]
  if (nrow(daa) != nrow(kegg_abundance)) {
    stop("The ALDEx2 GLM result is incomplete.", call. = FALSE)
  }
  required_p_columns <- c("p_values", "p_adjust")
  missing_p_columns <- setdiff(required_p_columns, colnames(daa))
  if (length(missing_p_columns)) {
    stop(
      "The ALDEx2 GLM result is missing p-value column(s): ",
      paste(missing_p_columns, collapse = ", "),
      call. = FALSE
    )
  }

  # Preserve the package-returned adjustment for provenance. Figure 6 uses one
  # four-treatment omnibus GLM p-value per pathway and does not calculate
  # post-hoc pairwise tests. The author additionally defined a conservative
  # decision rule motivated by the six possible pairs among four treatments:
  # multiply each raw omnibus p-value by six and cap it at one. The legacy-
  # compatible p_adjust column remains the value used for pathway selection.
  daa$raw_omnibus_p <- as.numeric(daa$p_values)
  daa$package_p_adjust <- as.numeric(daa$p_adjust)
  daa$package_adj_method <- if ("adj_method" %in% colnames(daa)) {
    as.character(daa$adj_method)
  } else {
    NA_character_
  }
  daa$author_sixfold_omnibus_p <- pmin(1, 6 * daa$raw_omnibus_p)
  daa$p_adjust <- daa$author_sixfold_omnibus_p
  daa$adj_method <- "author-defined sixfold conservative omnibus rule"
  daa$inference_scope <- "four-treatment omnibus only; no post-hoc pairwise tests"
  daa
}

load_offline_kegg_pathway_annotations <- function() {
  require_figure_6_packages()
  reference_environment <- new.env(parent = emptyenv())
  reference_file <- system.file(
    "extdata", "KO_reference.RData", package = "ggpicrust2"
  )
  if (!nzchar(reference_file) || !file.exists(reference_file)) {
    stop("ggpicrust2 KO_reference.RData is unavailable.", call. = FALSE)
  }
  load(reference_file, envir = reference_environment)
  ko_reference <- get("KO_reference", envir = reference_environment)

  pathway_id <- sub(
    ".*\\[PATH:(ko[0-9]{5})\\]$", "\\1", ko_reference$Pathway
  )
  has_pathway <- grepl("\\[PATH:ko[0-9]{5}\\]$", ko_reference$Pathway)
  annotations <- unique(data.frame(
    feature = pathway_id[has_pathway],
    pathway_name = sub(" \\[PATH:ko[0-9]{5}\\]$", "", ko_reference$Pathway[has_pathway]),
    pathway_class = paste0(
      ko_reference$PathwayL1[has_pathway], ";  ",
      ko_reference$PathwayL2[has_pathway]
    ),
    stringsAsFactors = FALSE
  ))
  duplicated_features <- unique(annotations$feature[duplicated(annotations$feature)])
  if (length(duplicated_features) > 0L) {
    # A pathway may occur repeatedly through its member KOs, but its pathway
    # name and hierarchy must be unique for deterministic local annotation.
    split_annotations <- split(annotations, annotations$feature)
    inconsistent <- names(Filter(function(x) {
      length(unique(x$pathway_name)) != 1L || length(unique(x$pathway_class)) != 1L
    }, split_annotations))
    if (length(inconsistent) > 0L) {
      stop(
        "Inconsistent offline KEGG annotations for: ",
        paste(inconsistent, collapse = ", "),
        call. = FALSE
      )
    }
    annotations <- do.call(rbind, lapply(split_annotations, function(x) x[1L, ]))
    rownames(annotations) <- NULL
  }
  annotations
}

annotate_figure_6_daa <- function(daa) {
  annotations <- load_offline_kegg_pathway_annotations()
  annotated <- merge(daa, annotations, by = "feature", all.x = TRUE, sort = FALSE)
  annotated <- annotated[match(daa$feature, annotated$feature), , drop = FALSE]
  annotated
}

prepare_figure_6_analysis <- function(
    project_root,
    batch_id,
    generation,
    seed = 20260809L) {
  input <- read_figure_6_input(project_root, batch_id, generation)
  abundance <- aggregate_ko_to_kegg(input)
  daa <- run_figure_6_daa(abundance, input$metadata, seed = seed)
  daa <- annotate_figure_6_daa(daa)
  list(
    input = input,
    abundance = abundance,
    daa = daa,
    seed = seed
  )
}

summarize_pathway_relative_abundance <- function(
    abundance,
    metadata,
    features) {
  missing_features <- setdiff(features, rownames(abundance))
  if (length(missing_features) > 0L) {
    stop(
      "Selected pathways missing from aggregated abundance: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
  }
  metadata <- metadata[match(colnames(abundance), metadata$sample_id), , drop = FALSE]
  if (anyNA(metadata$sample_id)) {
    stop("Could not align relative abundance values and metadata.", call. = FALSE)
  }

  values <- as.matrix(abundance)
  storage.mode(values) <- "double"
  sample_totals <- colSums(values)
  if (any(!is.finite(sample_totals)) || any(sample_totals <= 0)) {
    stop("One or more samples have invalid total KEGG abundance.", call. = FALSE)
  }
  relative <- sweep(values, 2L, sample_totals, "/")

  rows <- lapply(features, function(feature) {
    group_rows <- lapply(figure_6_treatment_levels, function(treatment) {
      current <- relative[
        feature,
        metadata$treatment == treatment,
        drop = TRUE
      ]
      data.frame(
        feature = feature,
        group = treatment,
        mean = mean(current),
        sd = stats::sd(current),
        n = length(current),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, group_rows)
  })
  summary <- do.call(rbind, rows)
  summary$group <- factor(summary$group, levels = figure_6_treatment_levels)
  rownames(summary) <- NULL
  summary
}

format_figure_6_p_value <- function(x) {
  formatC(x, format = "f", digits = 4L)
}

# Corrected multi-group replacement for the legacy pathway_errorbar3().
# It keeps the public arguments relevant to Figure 6, but aligns every table by
# explicit feature/group keys and never computes an undefined two-group fold
# change when four treatment groups are present.
pathway_errorbar3 <- function(
    abundance,
    daa_results_df,
    Group,
    ko_to_kegg = TRUE,
    p_values_threshold = 0.05,
    order = "pathway_class",
    select = NULL,
    p_value_bar = FALSE,
    colors = figure_6_treatment_colors,
    x_lab = "pathway_name",
    panel_letter = NULL,
    title = NULL,
    show_legend = TRUE) {
  require_figure_6_packages()
  if (!isTRUE(ko_to_kegg)) {
    stop("Figure 6 pathway_errorbar3() supports KO-to-KEGG data only.", call. = FALSE)
  }
  if (isTRUE(p_value_bar)) {
    warning(
      "p_value_bar is disabled for the four-group Figure 6 comparison.",
      call. = FALSE
    )
  }
  required <- c("feature", "p_adjust", "pathway_class", x_lab)
  if (!all(required %in% names(daa_results_df))) {
    stop("The annotated DAA table is missing plotting columns.", call. = FALSE)
  }
  if (length(Group) != ncol(abundance)) {
    stop("Group length must equal the number of abundance columns.", call. = FALSE)
  }
  Group <- factor(Group, levels = figure_6_treatment_levels)
  if (anyNA(Group)) stop("Group contains an invalid treatment label.", call. = FALSE)

  threshold_selected <- daa_results_df[
    is.finite(daa_results_df$p_adjust) &
      daa_results_df$p_adjust <= p_values_threshold,
    ,
    drop = FALSE
  ]
  if (!is.null(select)) {
    missing_select <- setdiff(select, daa_results_df$feature)
    if (length(missing_select) > 0L) {
      stop(
        "Requested pathways are absent from the DAA result: ",
        paste(missing_select, collapse = ", "),
        call. = FALSE
      )
    }
    selected <- daa_results_df[match(select, daa_results_df$feature), , drop = FALSE]
  } else {
    selected <- threshold_selected
  }
  if (nrow(selected) == 0L) {
    stop("No pathways pass the panel selection rule.", call. = FALSE)
  }
  if (anyNA(selected[[x_lab]]) || anyNA(selected$pathway_class)) {
    missing <- selected$feature[is.na(selected[[x_lab]]) | is.na(selected$pathway_class)]
    stop(
      "Missing offline KEGG annotation for selected pathway(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (identical(order, "pathway_class")) {
    selected <- selected[order(selected$pathway_class, selected$p_adjust), , drop = FALSE]
  } else if (identical(order, "p_values")) {
    selected <- selected[order(selected$p_adjust), , drop = FALSE]
  } else if (identical(order, "name")) {
    selected <- selected[order(selected[[x_lab]]), , drop = FALSE]
  } else if (is.numeric(order)) {
    selected <- selected[order, , drop = FALSE]
  }

  metadata <- data.frame(
    sample_id = colnames(abundance),
    treatment = Group,
    stringsAsFactors = FALSE
  )
  summary <- summarize_pathway_relative_abundance(
    abundance,
    metadata,
    selected$feature
  )
  summary <- merge(
    summary,
    selected[, c("feature", x_lab, "pathway_class", "p_adjust"), drop = FALSE],
    by = "feature",
    all.x = TRUE,
    sort = FALSE
  )
  summary <- summary[
    order(match(summary$feature, selected$feature), match(summary$group, figure_6_treatment_levels)),
    ,
    drop = FALSE
  ]

  # Reverse the feature levels so the first ordered pathway appears at top.
  feature_levels <- rev(selected$feature)
  summary$feature_factor <- factor(summary$feature, levels = feature_levels)
  selected$feature_factor <- factor(selected$feature, levels = feature_levels)
  class_levels <- unique(selected$pathway_class)
  class_palette <- c(
    "Genetic Information Processing;  Folding, sorting and degradation" = "#8A9FD1",
    "Genetic Information Processing;  Replication and repair" = "#FEE500",
    "Metabolism;  Amino acid metabolism" = "#F47D2B",
    "Metabolism;  Energy metabolism" = "#89288F",
    "Metabolism;  Metabolism of cofactors and vitamins" = "#208A42",
    "Metabolism;  Metabolism of other amino acids" = "#272E6A",
    "Metabolism;  Xenobiotics biodegradation and metabolism" = "#D51F26",
    "Metabolism;  Biosynthesis of other secondary metabolites" = "#D51F26"
  )
  class_colors <- unname(class_palette[class_levels])
  if (anyNA(class_colors)) {
    fallback <- c("#C06CAB", "#90D5E4", "#89C75F", "#F37B7D")
    class_colors[is.na(class_colors)] <- fallback[seq_len(sum(is.na(class_colors)))]
  }
  class_colors <- setNames(class_colors, class_levels)
  class_rows <- split(selected$feature, selected$pathway_class)
  class_background <- do.call(rbind, lapply(class_levels, function(pathway_class) {
    positions <- match(class_rows[[pathway_class]], feature_levels)
    data.frame(
      pathway_class = pathway_class,
      ymin = min(positions) - 0.5,
      ymax = max(positions) + 0.5,
      stringsAsFactors = FALSE
    )
  }))
  class_background$fill <- unname(class_colors[class_background$pathway_class])
  class_background$y <- (class_background$ymin + class_background$ymax) / 2
  class_background$text_size <- ifelse(
    class_background$pathway_class ==
      "Metabolism;  Biosynthesis of other secondary metabolites",
    2.60,
    2.85
  )
  class_background$text_x <- ifelse(
    class_background$pathway_class ==
      "Metabolism;  Biosynthesis of other secondary metabolites",
    0.30,
    0.35
  )

  common_theme <- ggplot2::theme_minimal(base_family = "sans", base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "#F1F1F1", color = NA)
    )

  selected$feature_y <- match(selected$feature, feature_levels)
  class_plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = class_background,
      ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = class_background$fill,
      alpha = 0.23,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = class_background,
      ggplot2::aes(
        x = text_x,
        y = y,
        label = pathway_class,
        size = text_size
      ),
      fontface = "bold",
      lineheight = 0.95
    ) +
    ggplot2::geom_text(
      data = selected,
      ggplot2::aes(x = 0.99, y = feature_y, label = .data[[x_lab]]),
      hjust = 1,
      fontface = "bold",
      size = 3.05
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_size_identity() +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0.5, length(feature_levels) + 0.5),
      clip = "off"
    ) +
    common_theme +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 3, 0)
    )

  x_labels <- if (max(summary$mean + summary$sd, na.rm = TRUE) < 1e-4) {
    function(x) formatC(x, format = "e", digits = 1L)
  } else {
    scales::label_number(accuracy = NULL)
  }

  bar_plot <- ggplot2::ggplot(
    summary,
    ggplot2::aes(x = mean, y = feature_factor, fill = group)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = 0, xmax = mean + sd),
      position = ggplot2::position_dodge(width = 0.78),
      width = 0.48,
      linewidth = 0.45,
      color = "black"
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.78),
      width = 0.78
    ) +
    ggplot2::scale_fill_manual(
      values = colors,
      limits = figure_6_treatment_levels,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::scale_y_discrete(
      labels = NULL,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.06)),
      labels = x_labels
    ) +
    ggplot2::labs(x = "Relative Abundance", y = NULL) +
    common_theme +
    ggplot2::theme(
      axis.line.x = ggplot2::element_line(linewidth = 0.45, color = "black"),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.45, color = "black"),
      axis.text.x = ggplot2::element_text(color = "black", face = "bold", size = 9),
      axis.text.y = ggplot2::element_text(color = "black", face = "bold", size = 9),
      axis.title.x = ggplot2::element_text(face = "bold", size = 10),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 2, 3, 0)
    )

  p_plot <- ggplot2::ggplot(
    selected,
    ggplot2::aes(x = 0.5, y = feature_factor, label = format_figure_6_p_value(p_adjust))
  ) +
    ggplot2::geom_text(fontface = "bold", size = 3.3) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(x = NULL, y = "adj. p-value") +
    common_theme +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(face = "bold", size = 10),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 1, 3, 0)
    )

  body_widths <- if (
    length(feature_levels) == 1L &&
      identical(selected$pathway_class[[1L]],
                "Genetic Information Processing;  Folding, sorting and degradation")
  ) c(1.00, 1.00, 0.16) else c(1.22, 1.00, 0.16)
  body <- class_plot + bar_plot + p_plot +
    patchwork::plot_layout(widths = body_widths)

  header <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0.5, y = 0.52, label = title,
      fontface = "bold", size = 8.0
    ) +
    ggplot2::annotate(
      "text", x = 0.002, y = 0.52, label = panel_letter,
      hjust = 0, fontface = "bold", size = 9.0
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  if (isTRUE(show_legend)) {
    legend_data <- data.frame(
      group = factor(figure_6_treatment_levels, levels = figure_6_treatment_levels),
      x = 0.875,
      y = c(0.82, 0.61, 0.40, 0.19),
      stringsAsFactors = FALSE
    )
    header <- header +
      ggplot2::geom_tile(
        data = legend_data,
        ggplot2::aes(x = x, y = y, fill = group),
        width = 0.040,
        height = 0.16
      ) +
      ggplot2::geom_text(
        data = legend_data,
        ggplot2::aes(x = x + 0.030, y = y, label = group),
        hjust = 0,
        fontface = "bold",
        size = 4.0
      ) +
      ggplot2::scale_fill_manual(
        values = colors,
        limits = figure_6_treatment_levels,
        drop = FALSE
      ) +
      ggplot2::theme(legend.position = "none")
  }
  header <- patchwork::wrap_elements(full = header)
  header_weight <- if (length(feature_levels) > 2L) 0.17 else 0.55
  plot <- header / body + patchwork::plot_layout(heights = c(header_weight, 1))

  attr(plot, "selected_daa") <- selected
  attr(plot, "relative_abundance_summary") <- summary
  plot
}

save_figure_6_panel <- function(plot, filename, width_px, height_px, dpi = 180) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width_px / dpi,
    height = height_px / dpi,
    units = "in",
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  invisible(filename)
}
