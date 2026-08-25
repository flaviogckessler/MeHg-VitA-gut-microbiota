# Reproducible refactor of
# scripts/legacy_scripts/heatmap_zscore_legacy.R.
# Panel A: significant taxon relative abundances in dams (F0).
#
# The legacy script is retained for provenance and is not a runtime dependency.
# This standalone copy uses canonical files organized under New_repo and
# reproduces the legacy logic:
# Kruskal-Wallis screening (p <= 0.05), treatment means, and row-wise z-scores.

# Load Packages ----
suppressPackageStartupMessages({
  library(graphics)
  library(grDevices)
  library(stats)
})

# - Repository paths
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata"),
      file.path(path, "results", "taxon_relative_abundance"),
      file.path(path, "figures", "Figure_4")
    )
    if (all(dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the New_repo project root.", call. = FALSE)
    }
    path <- parent
  }
}

treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
rank_levels <- c("p", "c", "o", "f", "g")

# - Canonical taxon-abundance input
read_taxon_abundance <- function(project_root, batch_id, filename) {
  abundance_file <- file.path(
    project_root, "results", "taxon_relative_abundance", filename
  )
  metadata_file <- file.path(
    project_root, "data", "metadata", batch_id, "sample_manifest.csv"
  )

  abundance <- read.csv(
    abundance_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )
  metadata <- read.csv(
    metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )
  included <- tolower(as.character(metadata$include_primary_analysis)) %in%
    c("true", "t", "1", "yes")
  metadata <- metadata[included, , drop = FALSE]

  required_columns <- c("tax_rank", "tax_name", "regex_match")
  if (!all(required_columns %in% names(abundance))) {
    stop("The abundance table is missing required taxon columns.", call. = FALSE)
  }
  analyzed_labels <- abundance$regex_match[
    !is.na(abundance$regex_match) & abundance$tax_rank %in% rank_levels
  ]
  if (anyDuplicated(analyzed_labels)) {
    stop("The analyzed taxonomic ranks contain duplicated labels.", call. = FALSE)
  }

  list(
    abundance = abundance,
    metadata = metadata,
    source_files = c(abundance = abundance_file, metadata = metadata_file)
  )
}

# - Statistics copied from the legacy heatmap workflow
calculate_taxon_stratum <- function(input, generation, sex = NULL) {
  metadata <- input$metadata[input$metadata$generation == generation, , drop = FALSE]
  if (!is.null(sex)) {
    metadata <- metadata[metadata$sex == sex, , drop = FALSE]
  }
  metadata$treatment <- factor(metadata$treatment, levels = treatment_levels)

  sample_ids <- metadata$sample_id[metadata$sample_id %in% names(input$abundance)]
  if (length(sample_ids) != nrow(metadata)) {
    missing <- setdiff(metadata$sample_id, names(input$abundance))
    stop(
      paste("Included samples missing from abundance table:", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]

  keep_taxa <- !is.na(input$abundance$regex_match) &
    input$abundance$tax_rank %in% rank_levels
  taxa <- input$abundance[keep_taxa, , drop = FALSE]
  values <- as.matrix(taxa[, sample_ids, drop = FALSE])
  storage.mode(values) <- "double"
  nonzero <- rowSums(values, na.rm = TRUE) > 0
  taxa <- taxa[nonzero, , drop = FALSE]
  values <- values[nonzero, , drop = FALSE]

  p_values <- apply(values, 1, function(current_values) {
    tryCatch(
      suppressWarnings(
        kruskal.test(current_values, metadata$treatment)$p.value
      ),
      error = function(e) NA_real_
    )
  })
  significant <- is.finite(p_values) & p_values <= 0.05
  taxa <- taxa[significant, , drop = FALSE]
  values <- values[significant, , drop = FALSE]
  p_values <- p_values[significant]

  ordering <- order(
    match(taxa$tax_rank, rank_levels),
    taxa$regex_match,
    method = "radix"
  )
  taxa <- taxa[ordering, , drop = FALSE]
  values <- values[ordering, , drop = FALSE]
  p_values <- p_values[ordering]

  group_means <- vapply(treatment_levels, function(treatment) {
    rowMeans(
      values[, metadata$treatment == treatment, drop = FALSE],
      na.rm = TRUE
    )
  }, numeric(nrow(values)))
  colnames(group_means) <- treatment_levels
  z_scores <- t(apply(group_means, 1, function(x) (x - mean(x)) / stats::sd(x)))
  colnames(z_scores) <- treatment_levels

  output <- data.frame(
    taxon = taxa$regex_match,
    rank = taxa$tax_rank,
    p_value = unname(p_values),
    stratum = if (is.null(sex)) "Dams" else sex,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (treatment in treatment_levels) {
    output[[paste0("mean_", treatment)]] <- group_means[, treatment]
    output[[paste0("z_", treatment)]] <- z_scores[, treatment]
  }

  attr(output, "sample_ids") <- sample_ids
  output
}

prepare_figure_4_a_data <- function(project_root = find_project_root()) {
  input <- read_taxon_abundance(
    project_root,
    "batch01_iontorrent_dams-f0_f1-pnd21",
    "batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv"
  )
  data <- calculate_taxon_stratum(input, generation = "F0")
  expected_taxa <- c(
    "g__[Eubacterium] ventriosum group",
    "g__[Ruminococcus] gauvreauii group",
    "g__Oribacterium"
  )
  if (!setequal(data$taxon, expected_taxa)) {
    stop("Panel A significant taxa differ from the final-figure specification.", call. = FALSE)
  }
  data <- data[match(expected_taxa, data$taxon), , drop = FALSE]
  attr(data, "source_files") <- input$source_files
  data
}

# - Plotting helpers
format_p_value <- function(p) {
  paste0("p = ", format(signif(p, 2), scientific = FALSE, trim = TRUE))
}

draw_figure_4_heatmap <- function(
    data,
    filename,
    panel_letter,
    title,
    width_px,
    height_px,
    margins_in,
    show_sex_labels = FALSE) {
  z_columns <- paste0("z_", treatment_levels)
  z_scores <- as.matrix(data[, z_columns, drop = FALSE])
  storage.mode(z_scores) <- "double"
  row_count <- nrow(z_scores)
  if (row_count == 0) stop("No significant taxa were found.", call. = FALSE)

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  png_type <- if (capabilities("cairo")) "cairo" else "windows"
  grDevices::png(
    filename,
    width = width_px,
    height = height_px,
    units = "px",
    res = 200,
    pointsize = 10,
    bg = "white",
    type = png_type
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(
    mai = margins_in,
    family = "sans",
    xaxs = "i",
    yaxs = "i"
  )
  palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(15)
  breaks <- seq(-1.5, 1.5, length.out = length(palette) + 1L)
  graphics::image(
    x = seq_along(treatment_levels),
    y = seq_len(row_count),
    z = t(z_scores[row_count:1, , drop = FALSE]),
    xlim = c(0.5, 4.5),
    ylim = c(0.5, row_count + 0.5),
    col = palette,
    breaks = breaks,
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = TRUE
  )

  row_y <- rev(seq_len(row_count))
  p_labels <- if ("p_label" %in% names(data)) {
    data$p_label
  } else {
    vapply(data$p_value, format_p_value, character(1))
  }
  graphics::text(
    x = 0.36,
    y = row_y,
    labels = p_labels,
    adj = 1,
    cex = 0.80,
    xpd = NA
  )
  graphics::text(
    x = 4.64,
    y = row_y,
    labels = data$taxon,
    adj = 0,
    cex = 0.94,
    xpd = NA
  )
  graphics::text(
    x = seq_along(treatment_levels),
    y = 0.27,
    labels = treatment_levels,
    srt = 90,
    adj = 1,
    cex = 0.94,
    xpd = NA
  )
  graphics::text(
    x = -0.82,
    y = row_count + 0.72,
    labels = "Kruskal-Wallis:",
    adj = 0,
    font = 2,
    cex = 0.82,
    xpd = NA
  )
  graphics::text(
    x = 4.18,
    y = row_count + 1.02,
    labels = title,
    font = 2,
    cex = 1.42,
    xpd = NA
  )
  graphics::text(
    x = graphics::grconvertX(0.02, from = "ndc", to = "user"),
    y = row_count + 1.35,
    labels = panel_letter,
    adj = 0,
    font = 2,
    cex = 1.52,
    xpd = NA
  )

  if (show_sex_labels && length(unique(data$stratum)) > 1L) {
    strata <- unique(data$stratum)
    for (stratum in strata) {
      index <- which(data$stratum == stratum)
      graphics::text(
        x = -1.22,
        y = mean(row_y[index]),
        labels = paste0(stratum, "s"),
        srt = 90,
        font = 2,
        cex = 1.18,
        xpd = NA
      )
    }
    boundaries <- cumsum(table(factor(data$stratum, levels = strata)))
    boundaries <- boundaries[-length(boundaries)]
    for (boundary in boundaries) {
      separator_y <- row_count - boundary + 0.5
      graphics::rect(
        xleft = 0.5,
        ybottom = separator_y - 0.10,
        xright = 4.5,
        ytop = separator_y + 0.10,
        col = "white",
        border = NA
      )
      graphics::segments(
        x0 = -1.48,
        y0 = separator_y,
        x1 = 10.8,
        y1 = separator_y,
        col = "#666666",
        lty = 2,
        lwd = 0.7,
        xpd = NA
      )
    }
  }

  invisible(filename)
}

# - Figure generation
generate_figure_4_a <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_4", "output")) {
  project_root <- find_project_root()
  data <- prepare_figure_4_a_data(project_root)
  file <- file.path(output_dir, "Figure4A_generated.png")
  draw_figure_4_heatmap(
    data = data,
    filename = file,
    panel_letter = "A",
    title = "Dams (F0)",
    width_px = 1000,
    height_px = 360,
    margins_in = c(0.965, 0.88, 0.34, 3.145),
    show_sex_labels = FALSE
  )
  invisible(list(file = file, data = data))
}

if (sys.nframe() == 0) {
  result <- generate_figure_4_a()
  message("Generated Figure 4 panel A: ", normalizePath(result$file, winslash = "/"))
}
