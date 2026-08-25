# Reproducible Figure S7 panels A-D.
#
# Figure S7 shows PCA of the untransformed PICRUSt KO prediction scores for
# Dams (F0), F1-PND21, F1-PND365, and F2-PND365. The supplied
# Figure_S7_FinalFigure.png was assembled manually; this script generates the
# four individual PNG panels used in that layout. The historical workflow is
# preserved in scripts/legacy_scripts/picrust_legacy.R and is not
# sourced at runtime.

# Load Packages ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(stats)
  library(utils)
})

figure_s7_required_packages <- "ggplot2"
figure_s7_missing_packages <- figure_s7_required_packages[
  !vapply(
    figure_s7_required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(figure_s7_missing_packages) > 0L) {
  stop(
    "Install the package required by Figure S7: ",
    paste(figure_s7_missing_packages, collapse = ", "),
    call. = FALSE
  )
}

figure_s7_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
figure_s7_sex_levels <- c("Female", "Male")
figure_s7_colors <- c(
  CN = "#F8766D",
  VitA = "#7CAE00",
  MeHg = "#00BFC4",
  `VitA+MeHg` = "#C77CFF"
)
figure_s7_fills <- paste0(figure_s7_colors, "33")
figure_s7_shapes <- c(Female = 16, Male = 17)

figure_s7_batch_dimensions <- list(
  `batch01_iontorrent_dams-f0_f1-pnd21` = c(6710L, 55L),
  `batch02_illumina-miseq_f1-f2_pnd365` = c(6026L, 50L)
)

figure_s7_panel_configurations <- list(
  A = list(
    label = "Dams (F0)",
    batch_id = "batch01_iontorrent_dams-f0_f1-pnd21",
    generation = "F0",
    collection_day_pnd = NA_integer_,
    expected_sample_count = 18L,
    expected_sexes = "Female",
    expected_variance = c(75.36748, 10.42092),
    orientation_samples = c(PC1 = "M9F", PC2 = "M9F"),
    orientation_signs = c(PC1 = 1, PC2 = 1)
  ),
  B = list(
    label = "F1-PND21",
    batch_id = "batch01_iontorrent_dams-f0_f1-pnd21",
    generation = "F1",
    collection_day_pnd = 21L,
    expected_sample_count = 34L,
    expected_sexes = figure_s7_sex_levels,
    expected_variance = c(58.63000, 21.72620),
    orientation_samples = c(PC1 = "42F", PC2 = "11F"),
    orientation_signs = c(PC1 = -1, PC2 = 1)
  ),
  C = list(
    label = "F1-PND365",
    batch_id = "batch02_illumina-miseq_f1-f2_pnd365",
    generation = "F1",
    collection_day_pnd = 365L,
    expected_sample_count = 24L,
    expected_sexes = figure_s7_sex_levels,
    expected_variance = c(87.06747, 8.16765),
    orientation_samples = c(PC1 = "t8", PC2 = "t5"),
    orientation_signs = c(PC1 = 1, PC2 = -1)
  ),
  D = list(
    label = "F2-PND365",
    batch_id = "batch02_illumina-miseq_f1-f2_pnd365",
    generation = "F2",
    collection_day_pnd = 365L,
    expected_sample_count = 24L,
    expected_sexes = figure_s7_sex_levels,
    expected_variance = c(87.59635, 7.66853),
    orientation_samples = c(PC1 = "t38", PC2 = "t47"),
    orientation_signs = c(PC1 = 1, PC2 = 1)
  )
)

find_figure_s7_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata", "sample_manifest_all_batches.csv"),
      file.path(path, "results", "functional_prediction"),
      file.path(
        path, "supplementary_figures", "Figure_S7", "output",
        "Figure_S7_FinalFigure.png"
      )
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the repository root for Figure S7.", call. = FALSE)
    }
    path <- parent
  }
}

figure_s7_batch_files <- function(project_root, batch_id) {
  list(
    ko_abundance = file.path(
      project_root,
      "results",
      "functional_prediction",
      batch_id,
      "kegg_ortholog_abundance.tsv"
    ),
    metadata = file.path(
      project_root,
      "data",
      "metadata",
      batch_id,
      "sample_manifest.csv"
    )
  )
}

read_figure_s7_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)

  signature <- readBin(connection, what = "raw", n = 8L)
  expected_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (!identical(signature, expected_signature)) {
    stop("Not a valid PNG file: ", filename, call. = FALSE)
  }

  chunk_length <- readBin(
    connection, what = "integer", n = 1L, size = 4L, endian = "big"
  )
  chunk_type <- rawToChar(readBin(connection, what = "raw", n = 4L))
  if (!identical(chunk_length, 13L) || !identical(chunk_type, "IHDR")) {
    stop("PNG does not begin with an IHDR chunk: ", filename, call. = FALSE)
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

figure_s7_primary_flag <- function(values, filename) {
  normalized <- tolower(trimws(as.character(values)))
  if (!all(normalized %in% c("true", "false"))) {
    stop(
      "Invalid include_primary_analysis value in: ",
      filename,
      call. = FALSE
    )
  }
  normalized == "true"
}

read_figure_s7_batch <- function(
    batch_id,
    project_root = find_figure_s7_project_root()) {
  files <- figure_s7_batch_files(project_root, batch_id)
  missing_files <- unlist(files, use.names = FALSE)[
    !file.exists(unlist(files, use.names = FALSE))
  ]
  if (length(missing_files) > 0L) {
    stop(
      "Missing Figure S7 source file(s): ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  ko_abundance <- utils::read.delim(
    files$ko_abundance,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fileEncoding = "UTF-8-BOM"
  )
  expected_dimensions <- figure_s7_batch_dimensions[[batch_id]]
  if (is.null(expected_dimensions) ||
      !identical(dim(ko_abundance), expected_dimensions)) {
    stop("Unexpected PICRUSt KO matrix dimensions for ", batch_id, ".", call. = FALSE)
  }
  if (!identical(names(ko_abundance)[1L], "ID")) {
    stop("The Figure S7 KO source must begin with an ID column.", call. = FALSE)
  }
  ko_ids <- sub("^ko:", "", trimws(ko_abundance$ID))
  if (any(!grepl("^K[0-9]{5}$", ko_ids)) || anyDuplicated(ko_ids)) {
    stop("The Figure S7 KO source contains invalid or duplicated IDs.", call. = FALSE)
  }

  metadata <- utils::read.csv(
    files$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  required_metadata <- c(
    "sample_id", "batch_id", "generation", "collection_day_pnd",
    "treatment", "sex", "include_primary_analysis"
  )
  if (!all(required_metadata %in% names(metadata))) {
    stop("The Figure S7 manifest is missing required columns.", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id) || anyNA(metadata$sample_id)) {
    stop("The Figure S7 manifest contains invalid sample IDs.", call. = FALSE)
  }
  if (!all(metadata$batch_id == batch_id)) {
    stop("The Figure S7 manifest contains an unexpected batch ID.", call. = FALSE)
  }
  if (!all(metadata$treatment %in% figure_s7_treatment_levels) ||
      !all(metadata$sex %in% figure_s7_sex_levels)) {
    stop("The Figure S7 manifest contains an unexpected treatment or sex.", call. = FALSE)
  }
  metadata$include_primary_analysis <- figure_s7_primary_flag(
    metadata$include_primary_analysis,
    files$metadata
  )

  list(
    ko_abundance = ko_abundance,
    ko_ids = ko_ids,
    metadata = metadata,
    files = files
  )
}

prepare_figure_s7_panel_data <- function(
    panel,
    batch_sources = NULL,
    project_root = find_figure_s7_project_root()) {
  configuration <- figure_s7_panel_configurations[[panel]]
  if (is.null(configuration)) {
    stop("Unknown Figure S7 panel: ", panel, call. = FALSE)
  }

  if (is.null(batch_sources)) {
    batch_sources <- read_figure_s7_batch(
      configuration$batch_id,
      project_root
    )
  }
  metadata <- batch_sources$metadata
  metadata_day <- suppressWarnings(as.integer(metadata$collection_day_pnd))
  day_matches <- if (is.na(configuration$collection_day_pnd)) {
    is.na(metadata_day)
  } else {
    metadata_day == configuration$collection_day_pnd
  }
  selected_metadata <- metadata[
    metadata$include_primary_analysis &
      metadata$generation == configuration$generation &
      day_matches,
    ,
    drop = FALSE
  ]
  if (nrow(selected_metadata) != configuration$expected_sample_count) {
    stop(
      "Unexpected sample count for Figure S7 panel ", panel, ".",
      call. = FALSE
    )
  }
  if (!setequal(unique(selected_metadata$sex), configuration$expected_sexes)) {
    stop("Unexpected sex composition for Figure S7 panel ", panel, ".", call. = FALSE)
  }
  treatment_counts <- table(
    factor(selected_metadata$treatment, levels = figure_s7_treatment_levels)
  )
  if (any(treatment_counts < 3L)) {
    stop(
      "Each Figure S7 treatment ellipse requires at least three samples.",
      call. = FALSE
    )
  }

  sample_ids <- selected_metadata$sample_id
  missing_samples <- setdiff(sample_ids, names(batch_sources$ko_abundance))
  if (length(missing_samples) > 0L) {
    stop(
      "Figure S7 samples absent from the KO source: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }
  ko_matrix <- as.matrix(
    batch_sources$ko_abundance[, sample_ids, drop = FALSE]
  )
  storage.mode(ko_matrix) <- "double"
  rownames(ko_matrix) <- batch_sources$ko_ids
  if (anyNA(ko_matrix) || any(!is.finite(ko_matrix)) || any(ko_matrix < 0)) {
    stop("Figure S7 KO scores must be finite and non-negative.", call. = FALSE)
  }

  pca <- stats::prcomp(t(ko_matrix), center = TRUE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  if (length(variance) < 2L ||
      any(abs(variance[1:2] - configuration$expected_variance) > 1e-3)) {
    stop(
      "PCA variance differs from the Figure S7 reference for panel ",
      panel,
      ".",
      call. = FALSE
    )
  }

  scores <- pca$x[, 1:2, drop = FALSE]
  for (component in colnames(scores)) {
    anchor <- configuration$orientation_samples[[component]]
    desired_sign <- configuration$orientation_signs[[component]]
    if (!anchor %in% rownames(scores) || scores[anchor, component] == 0) {
      stop("Invalid PCA orientation anchor for panel ", panel, ".", call. = FALSE)
    }
    if (sign(scores[anchor, component]) != desired_sign) {
      scores[, component] <- -scores[, component]
      pca$rotation[, component] <- -pca$rotation[, component]
    }
  }
  if (any(abs(colMeans(scores)) > 1e-6 * max(abs(scores)))) {
    stop("Figure S7 PCA scores are not centered.", call. = FALSE)
  }

  selected_metadata <- selected_metadata[
    match(rownames(scores), selected_metadata$sample_id),
    ,
    drop = FALSE
  ]
  if (anyNA(selected_metadata$sample_id)) {
    stop("Figure S7 PCA scores could not be aligned to metadata.", call. = FALSE)
  }
  score_data <- data.frame(
    sample_id = selected_metadata$sample_id,
    treatment = factor(
      selected_metadata$treatment,
      levels = figure_s7_treatment_levels
    ),
    sex = factor(selected_metadata$sex, levels = figure_s7_sex_levels),
    PC1 = scores[, "PC1"],
    PC2 = scores[, "PC2"],
    stringsAsFactors = FALSE
  )

  ellipse_covariances <- lapply(
    split(score_data[, c("PC1", "PC2")], score_data$treatment),
    stats::cov
  )
  if (any(vapply(
    ellipse_covariances,
    function(x) any(!is.finite(x)) || det(x) <= 0,
    logical(1)
  ))) {
    stop("A Figure S7 treatment ellipse is singular.", call. = FALSE)
  }

  attr(score_data, "pca") <- pca
  attr(score_data, "variance") <- variance
  attr(score_data, "source_files") <- batch_sources$files
  score_data
}

make_figure_s7_panel <- function(data, panel) {
  configuration <- figure_s7_panel_configurations[[panel]]
  variance <- attr(data, "variance")

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = PC1, y = PC2, colour = treatment)
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(fill = treatment),
      type = "norm",
      level = 0.68,
      geom = "polygon",
      segments = 51L,
      linewidth = 1,
      show.legend = c(colour = TRUE, fill = FALSE, shape = FALSE)
    )

  if (identical(configuration$expected_sexes, "Female")) {
    plot <- plot + ggplot2::geom_point(size = 2.5)
  } else {
    plot <- plot + ggplot2::geom_point(
      ggplot2::aes(shape = sex),
      size = 2.5
    )
  }

  plot <- plot +
    ggplot2::scale_colour_manual(
      name = "Group",
      values = figure_s7_colors,
      breaks = figure_s7_treatment_levels,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = figure_s7_fills,
      breaks = figure_s7_treatment_levels,
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", variance[[1L]]),
      y = sprintf("PC2 (%.1f%%)", variance[[2L]])
    ) +
    ggplot2::theme_bw() +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(fill = "#EBEBEB", shape = 16, size = 2.5)
      )
    ) +
    ggplot2::theme(
      legend.direction = "vertical",
      legend.position = "right",
      legend.key = ggplot2::element_rect(fill = "white", colour = NA)
    )

  if (!identical(configuration$expected_sexes, "Female")) {
    plot <- plot + ggplot2::scale_shape_manual(
      name = "Sex",
      values = figure_s7_shapes,
      breaks = figure_s7_sex_levels,
      drop = FALSE
    )
  }
  plot
}

generate_figure_s7_ad <- function(
    output_dir = file.path(
      find_figure_s7_project_root(),
      "supplementary_figures",
      "Figure_S7",
      "output"
    ),
    project_root = find_figure_s7_project_root()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  batch_ids <- unique(vapply(
    figure_s7_panel_configurations,
    `[[`,
    character(1),
    "batch_id"
  ))
  batch_sources <- setNames(
    lapply(batch_ids, read_figure_s7_batch, project_root = project_root),
    batch_ids
  )
  panel_data <- lapply(names(figure_s7_panel_configurations), function(panel) {
    configuration <- figure_s7_panel_configurations[[panel]]
    prepare_figure_s7_panel_data(
      panel,
      batch_sources = batch_sources[[configuration$batch_id]],
      project_root = project_root
    )
  })
  names(panel_data) <- names(figure_s7_panel_configurations)
  plots <- Map(make_figure_s7_panel, panel_data, names(panel_data))

  files <- setNames(
    file.path(
      output_dir,
      paste0("FigureS7", names(figure_s7_panel_configurations), "_generated.png")
    ),
    names(figure_s7_panel_configurations)
  )
  for (panel in names(files)) {
    ggplot2::ggsave(
      filename = files[[panel]],
      plot = plots[[panel]],
      device = "png",
      width = 6,
      height = 3,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }

  generated_dimensions <- vapply(
    files,
    read_figure_s7_png_dimensions,
    integer(2)
  )
  if (!all(generated_dimensions["width", ] == 1800L) ||
      !all(generated_dimensions["height", ] == 900L)) {
    stop("One or more Figure S7 outputs has unexpected dimensions.", call. = FALSE)
  }

  reference_file <- file.path(output_dir, "Figure_S7_FinalFigure.png")
  if (!file.exists(reference_file)) {
    stop("The supplied Figure S7 reference image is missing.", call. = FALSE)
  }
  reference_dimensions <- read_figure_s7_png_dimensions(reference_file)
  if (!identical(unname(reference_dimensions), c(2082L, 1179L))) {
    stop("The supplied Figure S7 reference dimensions changed.", call. = FALSE)
  }

  summary <- do.call(rbind, lapply(names(panel_data), function(panel) {
    variance <- attr(panel_data[[panel]], "variance")
    data.frame(
      panel = panel,
      cohort = figure_s7_panel_configurations[[panel]]$label,
      samples = nrow(panel_data[[panel]]),
      PC1_percent = variance[[1L]],
      PC2_percent = variance[[2L]],
      output = basename(files[[panel]]),
      stringsAsFactors = FALSE
    )
  }))

  invisible(list(
    files = files,
    plots = plots,
    data = panel_data,
    summary = summary,
    reference_file = reference_file,
    reference_dimensions = reference_dimensions,
    generated_dimensions = generated_dimensions
  ))
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s7_ad()
  message(
    "Generated and validated Figure S7 panels A-D in: ",
    normalizePath(dirname(result$files[[1L]]), winslash = "/")
  )
  print(result$summary, row.names = FALSE, digits = 7)
}
