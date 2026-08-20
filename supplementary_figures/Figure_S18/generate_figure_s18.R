# Reproducible Figure S18.
#
# This script generates one 4-row by 4-column figure.
# figure. Rows are K00520, K08363, K08364, and K08365; columns are F0-Dams,
# F1-PND21, F1-PND365, and F2-PND365. Sexes are pooled for treatment boxplots,
# while point shapes retain the sex information. Treatment effects are evaluated
# separately in every dataset-by-KO cell using one-way ANOVA followed by all
# six Tukey HSD pairwise comparisons. Only Tukey-adjusted p < 0.05 results are
# annotated.

# Load Packages ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(grDevices)
  library(grid)
  library(stats)
})

required_packages <- "ggplot2"
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the package required by Figure S18: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

find_figure_s18_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata", "sample_manifest_all_batches.csv"),
      file.path(path, "results", "functional_prediction"),
      file.path(path, "supplementary_figures", "Figure_S18")
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the repository root for Figure S18.", call. = FALSE)
    }
    path <- parent
  }
}

read_figure_s18_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)

  signature <- readBin(connection, what = "raw", n = 8L)
  expected_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (!identical(signature, expected_signature)) {
    stop("Figure S18 output is not a valid PNG file: ", filename, call. = FALSE)
  }

  chunk_length <- readBin(
    connection, what = "integer", n = 1L, size = 4L, endian = "big"
  )
  chunk_type <- rawToChar(readBin(connection, what = "raw", n = 4L))
  if (!identical(chunk_length, 13L) || !identical(chunk_type, "IHDR")) {
    stop("Figure S18 PNG does not begin with an IHDR chunk.", call. = FALSE)
  }

  c(
    width = readBin(
      connection, what = "integer", n = 1L, size = 4L, endian = "big"
    ),
    height = readBin(
      connection, what = "integer", n = 1L, size = 4L, endian = "big"
    )
  )
}

figure_s18_configurations <- list(
  F0_DAMS = list(
    label = "F0-Dams",
    batch_id = "batch01_iontorrent_dams-f0_f1-pnd21",
    generation = "F0",
    collection_day_pnd = NA_integer_,
    expected_sample_count = 18L,
    expected_matrix_dimensions = c(6710L, 55L),
    expected_sexes = "Female"
  ),
  F1_PND21 = list(
    label = "F1-PND21",
    batch_id = "batch01_iontorrent_dams-f0_f1-pnd21",
    generation = "F1",
    collection_day_pnd = 21L,
    expected_sample_count = 34L,
    expected_matrix_dimensions = c(6710L, 55L),
    expected_sexes = c("Female", "Male")
  ),
  F1_PND365 = list(
    label = "F1-PND365",
    batch_id = "batch02_illumina-miseq_f1-f2_pnd365",
    generation = "F1",
    collection_day_pnd = 365L,
    expected_sample_count = 24L,
    expected_matrix_dimensions = c(6026L, 50L),
    expected_sexes = c("Female", "Male")
  ),
  F2_PND365 = list(
    label = "F2-PND365",
    batch_id = "batch02_illumina-miseq_f1-f2_pnd365",
    generation = "F2",
    collection_day_pnd = 365L,
    expected_sample_count = 24L,
    expected_matrix_dimensions = c(6026L, 50L),
    expected_sexes = c("Female", "Male")
  )
)

figure_s18_ko_ids <- c(
  "K00520",
  "K08363",
  "K08364",
  "K08365"
)

figure_s18_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
figure_s18_sex_levels <- c("Female", "Male")

figure_s18_source_files <- function(project_root, configuration) {
  list(
    ko_abundance = file.path(
      project_root,
      "results", "functional_prediction", configuration$batch_id,
      "kegg_ortholog_abundance.tsv"
    ),
    metadata = file.path(
      project_root,
      "data", "metadata", configuration$batch_id,
      "sample_manifest.csv"
    ),
    kegg_query = file.path(
      project_root,
      "data", "reference", "kegg",
      "kegg_api_find_ko_mercuric_orthologs.tsv"
    )
  )
}

prepare_figure_s18_data <- function(
    cohort_id,
    project_root = find_figure_s18_project_root()) {
  configuration <- figure_s18_configurations[[cohort_id]]
  if (is.null(configuration)) {
    stop("Unknown Figure S18 cohort: ", cohort_id, call. = FALSE)
  }

  files <- figure_s18_source_files(project_root, configuration)
  missing_files <- unlist(files)[!file.exists(unlist(files))]
  if (length(missing_files) > 0L) {
    stop(
      "Missing Figure S18 source file(s): ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  kegg_query <- read.delim(
    files$kegg_query,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (!identical(names(kegg_query), c("KO", "Protein", "Description"))) {
    stop("Unexpected KEGG query schema for Figure S18.", call. = FALSE)
  }
  kegg_query$ko_id <- sub("^ko:", "", trimws(kegg_query$KO))
  kegg_query$Protein <- trimws(kegg_query$Protein)
  kegg_query$Description <- trimws(kegg_query$Description)
  requested_ids <- figure_s18_ko_ids
  if (!all(requested_ids %in% kegg_query$ko_id)) {
    stop("One or more Figure S18 KOs lack a KEGG annotation.", call. = FALSE)
  }

  metadata <- read.csv(
    files$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  required_metadata_columns <- c(
    "sample_id", "batch_id", "generation", "collection_day_pnd",
    "treatment", "sex", "include_primary_analysis"
  )
  if (!all(required_metadata_columns %in% names(metadata))) {
    stop("The Figure S18 metadata source is missing required columns.", call. = FALSE)
  }
  include_primary <- tolower(metadata$include_primary_analysis) == "true"
  metadata_day <- suppressWarnings(as.integer(metadata$collection_day_pnd))
  collection_day_matches <- if (is.na(configuration$collection_day_pnd)) {
    is.na(metadata_day)
  } else {
    metadata_day == configuration$collection_day_pnd
  }
  metadata <- metadata[
    include_primary &
      metadata$batch_id == configuration$batch_id &
      metadata$generation == configuration$generation &
      collection_day_matches,
    ,
    drop = FALSE
  ]
  if (nrow(metadata) != configuration$expected_sample_count) {
    stop(
      "Unexpected sample count for the ", configuration$label, " cohort.",
      call. = FALSE
    )
  }
  if (anyDuplicated(metadata$sample_id) || anyNA(metadata$sample_id)) {
    stop("Invalid sample IDs in the Figure S18 cohort.", call. = FALSE)
  }

  cohort_counts <- table(
    factor(metadata$sex, levels = configuration$expected_sexes),
    factor(metadata$treatment, levels = figure_s18_treatment_levels)
  )
  if (any(cohort_counts < 3L)) {
    stop("Every sex-by-treatment cell must contain at least 3 samples.", call. = FALSE)
  }

  ko_abundance <- read.delim(
    files$ko_abundance,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (!identical(names(ko_abundance)[1], "ID")) {
    stop("The Figure S18 KO source must begin with an ID column.", call. = FALSE)
  }
  if (!identical(dim(ko_abundance), configuration$expected_matrix_dimensions)) {
    stop("Unexpected PICRUSt2 matrix dimensions for Figure S18.", call. = FALSE)
  }
  ko_abundance$ko_id <- sub("^ko:", "", ko_abundance$ID)
  if (anyDuplicated(ko_abundance$ko_id)) {
    stop("Duplicated KO IDs in the Figure S18 source.", call. = FALSE)
  }
  if (!all(metadata$sample_id %in% names(ko_abundance))) {
    stop("One or more Figure S18 samples are absent from the KO matrix.", call. = FALSE)
  }
  if (!all(requested_ids %in% ko_abundance$ko_id)) {
    stop("Every Figure S18 KO must be present in both KO matrices.", call. = FALSE)
  }

  selected <- ko_abundance[
    match(requested_ids, ko_abundance$ko_id),
    metadata$sample_id,
    drop = FALSE
  ]
  abundance_matrix <- t(as.matrix(selected))
  storage.mode(abundance_matrix) <- "numeric"
  colnames(abundance_matrix) <- requested_ids
  if (anyNA(abundance_matrix) || any(!is.finite(abundance_matrix)) ||
      any(abundance_matrix < 0)) {
    stop("Figure S18 abundance profiles must be finite and non-negative.", call. = FALSE)
  }

  annotations <- kegg_query[
    match(requested_ids, kegg_query$ko_id),
    c("ko_id", "Protein", "Description"),
    drop = FALSE
  ]
  names(annotations) <- c("ko_id", "protein", "description")
  annotations$ko_label <- paste0(
    annotations$ko_id, " (", annotations$protein, ")"
  )

  panel_data <- do.call(
    rbind,
    lapply(requested_ids, function(ko_id) {
      data.frame(
        sample_id = metadata$sample_id,
        cohort_id = cohort_id,
        cohort_label = configuration$label,
        treatment = metadata$treatment,
        sex = metadata$sex,
        ko_id = ko_id,
        ko_label = annotations$ko_label[annotations$ko_id == ko_id],
        picrust_score = abundance_matrix[, ko_id],
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(panel_data) <- NULL
  panel_data$cohort_id <- factor(
    panel_data$cohort_id,
    levels = names(figure_s18_configurations)
  )
  panel_data$cohort_label <- factor(
    panel_data$cohort_label,
    levels = vapply(
      figure_s18_configurations, `[[`, character(1), "label"
    )
  )
  panel_data$treatment <- factor(
    panel_data$treatment,
    levels = figure_s18_treatment_levels
  )
  panel_data$sex <- factor(
    panel_data$sex,
    levels = figure_s18_sex_levels
  )
  panel_data$ko_id <- factor(panel_data$ko_id, levels = requested_ids)
  panel_data$ko_label <- factor(
    panel_data$ko_label,
    levels = annotations$ko_label
  )

  list(
    cohort_id = cohort_id,
    configuration = configuration,
    data = panel_data,
    annotations = annotations,
    cohort_counts = cohort_counts,
    files = files
  )
}

figure_s18_significance <- function(p_adjusted) {
  ifelse(
    p_adjusted < 0.0001, "****",
    ifelse(
      p_adjusted < 0.001, "***",
      ifelse(p_adjusted < 0.01, "**", ifelse(p_adjusted < 0.05, "*", "ns"))
    )
  )
}

compute_figure_s18_statistics <- function(plot_data) {
  cell_keys <- unique(plot_data[c("cohort_id", "cohort_label", "ko_id", "ko_label")])
  cell_keys <- cell_keys[order(cell_keys$ko_id, cell_keys$cohort_id), , drop = FALSE]

  omnibus_rows <- vector("list", nrow(cell_keys))
  posthoc_rows <- vector("list", nrow(cell_keys))
  for (index in seq_len(nrow(cell_keys))) {
    key <- cell_keys[index, , drop = FALSE]
    cell_data <- plot_data[
      plot_data$cohort_id == key$cohort_id & plot_data$ko_id == key$ko_id,
      ,
      drop = FALSE
    ]
    model <- stats::aov(picrust_score ~ treatment, data = cell_data)
    omnibus_p <- summary(model)[[1]][["Pr(>F)"]][1]
    tukey <- as.data.frame(stats::TukeyHSD(model, "treatment")$treatment)
    comparisons <- rownames(tukey)
    comparison_parts <- strsplit(comparisons, "-", fixed = TRUE)
    treatment_1 <- vapply(comparison_parts, `[[`, character(1), 2L)
    treatment_2 <- vapply(comparison_parts, `[[`, character(1), 1L)

    omnibus_rows[[index]] <- data.frame(
      cohort_id = as.character(key$cohort_id),
      cohort_label = as.character(key$cohort_label),
      ko_id = as.character(key$ko_id),
      ko_label = as.character(key$ko_label),
      anova_treatment_p = omnibus_p,
      stringsAsFactors = FALSE
    )
    posthoc_rows[[index]] <- data.frame(
      cohort_id = as.character(key$cohort_id),
      cohort_label = as.character(key$cohort_label),
      ko_id = as.character(key$ko_id),
      ko_label = as.character(key$ko_label),
      comparison = comparisons,
      treatment_1 = treatment_1,
      treatment_2 = treatment_2,
      difference = tukey$diff,
      confidence_low = tukey$lwr,
      confidence_high = tukey$upr,
      p_adjusted = tukey[["p adj"]],
      stringsAsFactors = FALSE
    )
  }

  omnibus <- do.call(rbind, omnibus_rows)
  posthoc <- do.call(rbind, posthoc_rows)
  posthoc$significance <- figure_s18_significance(posthoc$p_adjusted)
  significant <- posthoc[posthoc$significance != "ns", , drop = FALSE]

  if (nrow(significant) > 0L) {
    significant$x_start <- match(
      significant$treatment_1, figure_s18_treatment_levels
    )
    significant$x_end <- match(
      significant$treatment_2, figure_s18_treatment_levels
    )
    significant$x_left <- pmin(significant$x_start, significant$x_end)
    significant$x_right <- pmax(significant$x_start, significant$x_end)
    significant <- significant[
      order(
        significant$ko_id,
        significant$cohort_id,
        significant$x_right - significant$x_left,
        significant$x_left
      ),
      ,
      drop = FALSE
    ]
    significant$annotation_rank <- ave(
      seq_len(nrow(significant)),
      interaction(significant$cohort_id, significant$ko_id, drop = TRUE),
      FUN = seq_along
    )

    cell_max <- aggregate(
      picrust_score ~ cohort_id + ko_id,
      data = plot_data,
      FUN = max
    )
    cell_range <- aggregate(
      picrust_score ~ cohort_id + ko_id,
      data = plot_data,
      FUN = function(values) diff(range(values))
    )
    names(cell_max)[3] <- "cell_max"
    names(cell_range)[3] <- "cell_range"
    significant <- merge(
      significant,
      merge(cell_max, cell_range, by = c("cohort_id", "ko_id")),
      by = c("cohort_id", "ko_id"),
      all.x = TRUE,
      sort = FALSE
    )
    significant$cell_range[significant$cell_range <= 0] <- pmax(
      abs(significant$cell_max[significant$cell_range <= 0]) * 0.1,
      1
    )
    significant$y_position <- significant$cell_max +
      significant$cell_range * (0.10 + 0.13 * (significant$annotation_rank - 1L))
    significant$tick_height <- significant$cell_range * 0.035
    significant$label_y <- significant$y_position + significant$cell_range * 0.025
  }

  list(
    omnibus = omnibus,
    posthoc = posthoc,
    significant = significant
  )
}

make_figure_s18 <- function(plot_data, statistics) {
  treatment_colors <- c(
    CN = "#98F5FF",
    VitA = "#B3EE3A",
    MeHg = "#CD2626",
    `VitA+MeHg` = "#4D4D4D"
  )
  no_significant_note <- if (nrow(statistics$significant) == 0L) {
    "No Tukey-adjusted pairwise comparison reached p < 0.05"
  } else {
    paste0(
      nrow(statistics$significant),
      " Tukey-adjusted significant comparison(s) annotated"
    )
  }

  set.seed(1804L)
  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = treatment, y = picrust_score, fill = treatment)
  ) +
    ggplot2::geom_boxplot(
      width = 0.70,
      linewidth = 0.48,
      outlier.shape = NA,
      colour = "black"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = sex),
      position = ggplot2::position_jitter(width = 0.10, height = 0, seed = 1804L),
      colour = "black",
      fill = "white",
      alpha = 0.62,
      size = 1.45,
      stroke = 0.35
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(ko_label),
      cols = ggplot2::vars(cohort_label),
      scales = "free_y",
      switch = "y",
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(values = treatment_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(Female = 21, Male = 24),
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE,
      labels = c(
        CN = "CN", VitA = "VitA", MeHg = "MeHg", `VitA+MeHg` = "VitA+\nMeHg"
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.03, 0.28))
    ) +
    ggplot2::labs(
      title = "Figure S18 - Mercuric-related KEGG ortholog abundance",
      subtitle = no_significant_note,
      x = "Exposure group",
      y = "PICRUSt2 KO abundance",
      fill = "Exposure",
      shape = "Sex",
      caption = paste0(
        "Sexes pooled for one-way treatment ANOVA and Tukey HSD; ",
        "* p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001 (adjusted)."
      )
    ) +
    ggplot2::theme_gray(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        colour = "black", size = 8.5, angle = 25, hjust = 1, vjust = 1
      ),
      axis.text.y = ggplot2::element_text(colour = "black", size = 8.5),
      axis.title = ggplot2::element_text(colour = "black", size = 11),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(
        colour = "black", size = 9, face = "bold"
      ),
      legend.text = ggplot2::element_text(colour = "black", size = 9),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.48, "cm"),
      panel.background = ggplot2::element_rect(fill = "gray94", colour = NA),
      panel.border = ggplot2::element_rect(
        fill = NA, colour = "black", linewidth = 0.5
      ),
      panel.grid.major = ggplot2::element_line(colour = "white", linewidth = 0.45),
      panel.grid.minor = ggplot2::element_line(colour = "white", linewidth = 0.22),
      panel.spacing = grid::unit(0.12, "cm"),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(
        fill = "#FFE6E6", colour = "#B22222", linewidth = 0.55
      ),
      strip.text.x = ggplot2::element_text(
        colour = "black", size = 10.5, face = "bold"
      ),
      strip.text.y.left = ggplot2::element_text(
        colour = "black", size = 9, face = "bold", angle = 0
      ),
      plot.title = ggplot2::element_text(
        colour = "black", face = "bold", size = 16, hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "black", face = "bold", size = 10, hjust = 0.5
      ),
      plot.caption = ggplot2::element_text(
        colour = "black", size = 9, hjust = 0.5
      ),
      plot.margin = ggplot2::margin(t = 8, r = 10, b = 8, l = 8)
    )

  significant <- statistics$significant
  if (nrow(significant) > 0L) {
    significant$cohort_label <- factor(
      significant$cohort_label,
      levels = levels(plot_data$cohort_label)
    )
    significant$ko_label <- factor(
      significant$ko_label,
      levels = levels(plot_data$ko_label)
    )
    figure <- figure +
      ggplot2::geom_segment(
        data = significant,
        ggplot2::aes(
          x = x_left, xend = x_right,
          y = y_position, yend = y_position
        ),
        inherit.aes = FALSE,
        linewidth = 0.55
      ) +
      ggplot2::geom_segment(
        data = significant,
        ggplot2::aes(
          x = x_left, xend = x_left,
          y = y_position, yend = y_position - tick_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.55
      ) +
      ggplot2::geom_segment(
        data = significant,
        ggplot2::aes(
          x = x_right, xend = x_right,
          y = y_position, yend = y_position - tick_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.55
      ) +
      ggplot2::geom_text(
        data = significant,
        ggplot2::aes(
          x = (x_left + x_right) / 2,
          y = label_y,
          label = significance
        ),
        inherit.aes = FALSE,
        size = 4.2,
        fontface = "bold",
        vjust = 0
      )
  }

  figure
}

generate_figure_s18 <- function(
    output_dir = file.path(
      find_figure_s18_project_root(),
      "supplementary_figures", "Figure_S18", "output"
    ),
    project_root = find_figure_s18_project_root()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  prepared <- lapply(
    names(figure_s18_configurations),
    prepare_figure_s18_data,
    project_root = project_root
  )
  names(prepared) <- names(figure_s18_configurations)
  plot_data <- do.call(rbind, lapply(prepared, `[[`, "data"))
  rownames(plot_data) <- NULL

  expected_observations <- sum(vapply(
    figure_s18_configurations,
    `[[`,
    integer(1),
    "expected_sample_count"
  )) * length(figure_s18_ko_ids)
  if (nrow(plot_data) != expected_observations) {
    stop("Unexpected number of plotted observations in Figure S18.", call. = FALSE)
  }
  if ("K19058" %in% as.character(plot_data$ko_id)) {
    stop("K19058 must not appear in the four-KO Figure S18 workflow.", call. = FALSE)
  }

  statistics <- compute_figure_s18_statistics(plot_data)
  expected_posthoc_rows <-
    length(figure_s18_configurations) *
    length(figure_s18_ko_ids) *
    choose(length(figure_s18_treatment_levels), 2L)
  if (nrow(statistics$posthoc) != expected_posthoc_rows) {
    stop("Figure S18 must evaluate all six treatment pairs in every cell.", call. = FALSE)
  }

  figure <- make_figure_s18(plot_data, statistics)
  output_file <- file.path(
    output_dir,
    "FigureS18_generated.png"
  )
  grDevices::png(
    filename = output_file,
    width = 2700,
    height = 2400,
    units = "px",
    res = 180,
    bg = "white"
  )
  print(figure)
  grDevices::dev.off()

  dimensions <- read_figure_s18_png_dimensions(output_file)
  if (!identical(unname(dimensions), c(2700L, 2400L))) {
    stop("Unexpected Figure S18 PNG dimensions.", call. = FALSE)
  }

  invisible(list(
    file = output_file,
    dimensions = dimensions,
    data = plot_data,
    omnibus = statistics$omnibus,
    posthoc = statistics$posthoc,
    significant = statistics$significant,
    prepared = prepared,
    plot = figure
  ))
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s18()
  message("Generated official Figure S18 output: ", result$file)
  message(
    "Dimensions: ", result$dimensions[["width"]], " x ",
    result$dimensions[["height"]], " px"
  )
  message(
    "Post-hoc comparisons: ", nrow(result$posthoc),
    " total; ", nrow(result$significant),
    " Tukey-adjusted p < 0.05"
  )
  print(result$omnibus, row.names = FALSE)
  if (nrow(result$significant) > 0L) {
    print(
      result$significant[
        , c("cohort_id", "ko_id", "comparison", "p_adjusted", "significance")
      ],
      row.names = FALSE
    )
  }
}
