# Reproducible refactor of
# scripts/legacy_scripts/Correlacoes_Metadados.R.
# Figure 5 panel D: sex-stratified association between fecal total mercury
# (THg) and the relative abundance of g__Dorea in F1 offspring at PND21.
#
# The legacy script is retained for provenance, and the supplied TIFF was used
# as a visual migration reference. This generator uses canonical repository
# paths, applies the primary-analysis exclusions from the sample manifest, and
# has no non-base R dependencies or runtime dependency on the reference TIFF.

# - Repository paths
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata"),
      file.path(path, "results", "taxon_relative_abundance"),
      file.path(path, "figures", "Figure_5")
    )
    if (all(dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the New_repo project root.", call. = FALSE)
    }
    path <- parent
  }
}

batch_id <- "batch01_iontorrent_dams-f0_f1-pnd21"
sex_levels <- c("Female", "Male")
sex_colors <- c(Female = "#E64B35", Male = "#4DBBD5")

read_repository_csv <- function(filename, check.names = TRUE) {
  read.csv(
    filename,
    check.names = check.names,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )
}

check_required_columns <- function(data, columns, source_name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(
      paste(source_name, "is missing columns:", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

# - Canonical inputs
prepare_figure_5_d_data <- function(project_root = find_project_root()) {
  source_files <- c(
    thg_metadata = file.path(
      project_root, "data", "metadata", batch_id, "metadata_thg_hplc.csv"
    ),
    sample_manifest = file.path(
      project_root, "data", "metadata", batch_id, "sample_manifest.csv"
    ),
    taxon_relative_abundance = file.path(
      project_root,
      "results",
      "taxon_relative_abundance",
      "batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv"
    )
  )
  missing_files <- source_files[!file.exists(source_files)]
  if (length(missing_files) > 0L) {
    stop(
      paste("Figure 5D source files are missing:", paste(missing_files, collapse = ", ")),
      call. = FALSE
    )
  }

  thg_metadata <- read_repository_csv(source_files[["thg_metadata"]], check.names = FALSE)
  manifest <- read_repository_csv(source_files[["sample_manifest"]])
  abundance <- read_repository_csv(
    source_files[["taxon_relative_abundance"]],
    check.names = FALSE
  )

  check_required_columns(
    thg_metadata,
    c("SampleID", "Grupo", "Sexo", "Geracao", "Hg"),
    "THg/HPLC metadata"
  )
  check_required_columns(
    manifest,
    c(
      "sample_id", "generation", "life_stage", "collection_day_pnd",
      "treatment", "sex", "include_primary_analysis", "exclusion_reason"
    ),
    "Canonical sample manifest"
  )
  check_required_columns(
    abundance,
    c("tax_rank", "tax_name", "regex_match"),
    "Canonical taxon-relative-abundance table"
  )

  if (anyDuplicated(thg_metadata$SampleID)) {
    stop("THg/HPLC metadata contain duplicated sample identifiers.", call. = FALSE)
  }
  if (anyDuplicated(manifest$sample_id)) {
    stop("The canonical sample manifest contains duplicated sample identifiers.", call. = FALSE)
  }
  if (!all(thg_metadata$Geracao == "F1")) {
    stop("Figure 5D THg/HPLC metadata must contain only F1 offspring.", call. = FALSE)
  }

  manifest_index <- match(thg_metadata$SampleID, manifest$sample_id)
  if (anyNA(manifest_index)) {
    stop(
      paste(
        "THg/HPLC samples missing from the canonical manifest:",
        paste(thg_metadata$SampleID[is.na(manifest_index)], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  source_manifest <- manifest[manifest_index, , drop = FALSE]
  if (any(source_manifest$generation != "F1")) {
    stop("THg/HPLC samples do not map exclusively to canonical F1 records.", call. = FALSE)
  }
  if (any(thg_metadata$Sexo != source_manifest$sex)) {
    stop("Sex labels differ between THg/HPLC metadata and the manifest.", call. = FALSE)
  }
  if (any(thg_metadata$Grupo != source_manifest$treatment)) {
    stop("Treatment labels differ between THg/HPLC metadata and the manifest.", call. = FALSE)
  }

  include_primary <- tolower(as.character(source_manifest$include_primary_analysis)) %in%
    c("true", "t", "1", "yes")
  included_ids <- thg_metadata$SampleID[include_primary]
  excluded_ids <- thg_metadata$SampleID[!include_primary]
  if (!setequal(excluded_ids, c("10F", "22F"))) {
    stop(
      "The Figure 5D manifest exclusions differ from the documented legacy exclusions.",
      call. = FALSE
    )
  }

  dorea <- abundance[abundance$regex_match == "g__Dorea", , drop = FALSE]
  if (nrow(dorea) != 1L || dorea$tax_rank != "g" || dorea$tax_name != "Dorea") {
    stop("Expected exactly one genus-level g__Dorea row.", call. = FALSE)
  }
  missing_abundance <- setdiff(included_ids, names(dorea))
  if (length(missing_abundance) > 0L) {
    stop(
      paste(
        "Included Figure 5D samples missing from the abundance table:",
        paste(missing_abundance, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  thg_metadata <- thg_metadata[include_primary, , drop = FALSE]
  data <- data.frame(
    sample_id = thg_metadata$SampleID,
    sex = factor(thg_metadata$Sexo, levels = sex_levels),
    treatment = thg_metadata$Grupo,
    thg_ug_g = as.numeric(thg_metadata$Hg),
    dorea_relative_abundance = as.numeric(dorea[1, thg_metadata$SampleID]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (anyNA(data$sex) || anyNA(data$thg_ug_g) || anyNA(data$dorea_relative_abundance)) {
    stop("Figure 5D analysis data contain missing or invalid values.", call. = FALSE)
  }
  if (!identical(as.integer(table(data$sex)), c(14L, 14L))) {
    stop("Figure 5D must contain 14 female and 14 male primary-analysis samples.", call. = FALSE)
  }

  attr(data, "excluded_sample_ids") <- excluded_ids
  attr(data, "source_files") <- source_files
  data
}

# - Statistics reproduced from ggpubr::stat_cor(method = "spearman")
calculate_figure_5_d_statistics <- function(data) {
  statistics <- do.call(rbind, lapply(sex_levels, function(current_sex) {
    current <- data[data$sex == current_sex, , drop = FALSE]
    test <- suppressWarnings(stats::cor.test(
      current$thg_ug_g,
      current$dorea_relative_abundance,
      method = "spearman",
      exact = FALSE
    ))
    data.frame(
      sex = current_sex,
      n = nrow(current),
      rho = unname(test$estimate),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  }))
  rownames(statistics) <- NULL

  # Numeric checkpoints from the supplied legacy TIFF/workflow. Loose rounding
  # tolerances detect source drift while retaining the full recomputed values.
  expected <- data.frame(
    sex = sex_levels,
    n = c(14L, 14L),
    rho = c(0.81834, 0.65915),
    p_value = c(0.0003465, 0.01034),
    stringsAsFactors = FALSE
  )
  if (!identical(statistics$sex, expected$sex) ||
      !identical(statistics$n, expected$n) ||
      any(abs(statistics$rho - expected$rho) > 5e-5) ||
      any(abs(statistics$p_value - expected$p_value) > 5e-5)) {
    stop("Figure 5D Spearman statistics differ from the reference checkpoints.", call. = FALSE)
  }

  statistics
}

calculate_confidence_band <- function(data, grid_length = 200L) {
  fit <- stats::lm(dorea_relative_abundance ~ thg_ug_g, data = data)
  grid <- data.frame(
    thg_ug_g = seq(min(data$thg_ug_g), max(data$thg_ug_g), length.out = grid_length)
  )
  interval <- stats::predict(fit, newdata = grid, interval = "confidence", level = 0.95)
  data.frame(
    thg_ug_g = grid$thg_ug_g,
    fit = interval[, "fit"],
    lower = interval[, "lwr"],
    upper = interval[, "upr"]
  )
}

format_p_value <- function(p_value) {
  format(signif(p_value, 2), scientific = FALSE, trim = TRUE)
}

# - Reference-matched base-R plotting
draw_figure_5_d <- function(data, statistics, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  png_type <- if (capabilities("cairo")) "cairo" else "windows"
  grDevices::png(
    filename,
    width = 825,
    height = 600,
    units = "px",
    res = 150,
    pointsize = 12,
    bg = "white",
    type = png_type
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(
    mai = c(0.52, 0.77, 0.085, 1.18),
    mgp = c(1.75, 0.42, 0),
    tcl = -0.27,
    family = "sans",
    xaxs = "i",
    yaxs = "i"
  )
  x_limits <- c(-0.135, 2.75)
  y_limits <- c(-0.00255, 0.0133)
  graphics::plot.new()
  graphics::plot.window(xlim = x_limits, ylim = y_limits)

  bands <- lapply(sex_levels, function(current_sex) {
    calculate_confidence_band(data[data$sex == current_sex, , drop = FALSE])
  })
  names(bands) <- sex_levels

  # ggpubr adds the regression smoother after the scatter layer. Draw each
  # smoother group in factor order to preserve the same overlap appearance.
  for (current_sex in sex_levels) {
    current <- data[data$sex == current_sex, , drop = FALSE]
    graphics::points(
      current$thg_ug_g,
      current$dorea_relative_abundance,
      pch = 19,
      cex = 0.72,
      col = sex_colors[[current_sex]]
    )
  }
  for (current_sex in sex_levels) {
    band <- bands[[current_sex]]
    graphics::polygon(
      c(band$thg_ug_g, rev(band$thg_ug_g)),
      c(band$lower, rev(band$upper)),
      col = grDevices::adjustcolor(sex_colors[[current_sex]], alpha.f = 0.40),
      border = NA
    )
    graphics::lines(
      band$thg_ug_g,
      band$fit,
      col = sex_colors[[current_sex]],
      lwd = 1.45
    )
  }

  graphics::axis(
    1,
    at = 0:2,
    labels = 0:2,
    cex.axis = 1.00,
    lwd = 0,
    lwd.ticks = 1.05
  )
  y_ticks <- c(0, 0.004, 0.008, 0.012)
  graphics::axis(
    2,
    at = y_ticks,
    labels = sprintf("%.3f", y_ticks),
    las = 1,
    cex.axis = 1.00,
    lwd = 0,
    lwd.ticks = 1.05
  )
  graphics::segments(
    x_limits[1], y_limits[1], x_limits[2], y_limits[1],
    lwd = 1.05,
    xpd = FALSE
  )
  graphics::segments(
    x_limits[1], y_limits[1], x_limits[1], y_limits[2],
    lwd = 1.05,
    xpd = FALSE
  )
  graphics::mtext("THg, \u00b5g/g", side = 1, line = 1.55, cex = 1.08)
  graphics::mtext(
    expression(italic(Dorea~sp.)),
    side = 2,
    line = 2.40,
    cex = 1.08
  )

  label_positions <- data.frame(
    sex = sex_levels,
    x = c(0.15, 1.60),
    y = c(0.0042, 0.0075),
    stringsAsFactors = FALSE
  )
  for (current_sex in sex_levels) {
    current_statistics <- statistics[statistics$sex == current_sex, , drop = FALSE]
    position <- label_positions[label_positions$sex == current_sex, , drop = FALSE]
    correlation_label <- bquote(
      italic(R) == .(round(current_statistics$rho, 2)) * "," ~~
        italic(p) == .(as.numeric(format_p_value(current_statistics$p_value)))
    )
    graphics::text(
      position$x,
      position$y,
      labels = correlation_label,
      adj = c(0, 0.5),
      cex = 0.91,
      col = sex_colors[[current_sex]]
    )
  }

  # Draw a ggplot-style combined point/line/ribbon legend in the right margin.
  legend_x <- 2.96
  graphics::text(
    legend_x,
    0.00665,
    labels = "Sex",
    adj = c(0, 0.5),
    cex = 1.04,
    xpd = NA
  )
  legend_y <- c(Female = 0.00525, Male = 0.00416)
  for (current_sex in sex_levels) {
    current_y <- legend_y[[current_sex]]
    graphics::rect(
      legend_x,
      current_y - 0.00058,
      legend_x + 0.205,
      current_y + 0.00058,
      col = grDevices::adjustcolor(sex_colors[[current_sex]], alpha.f = 0.40),
      border = NA,
      xpd = NA
    )
    graphics::segments(
      legend_x + 0.035,
      current_y,
      legend_x + 0.170,
      current_y,
      col = sex_colors[[current_sex]],
      lwd = 1.45,
      xpd = NA
    )
    graphics::points(
      legend_x + 0.102,
      current_y,
      pch = 19,
      cex = 0.68,
      col = sex_colors[[current_sex]],
      xpd = NA
    )
    graphics::text(
      legend_x + 0.265,
      current_y,
      labels = current_sex,
      adj = c(0, 0.5),
      cex = 0.93,
      xpd = NA
    )
  }

  invisible(filename)
}

# - Output validation
read_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)
  header <- readBin(connection, what = "raw", n = 24L)
  signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(header) < 24L || !identical(header[1:8], signature)) {
    stop(paste("Invalid PNG signature:", filename), call. = FALSE)
  }
  unsigned_32 <- function(bytes) {
    sum(as.numeric(as.integer(bytes)) * 256^(3:0))
  }
  c(width = unsigned_32(header[17:20]), height = unsigned_32(header[21:24]))
}

validate_figure_5_d_output <- function(result, project_root = find_project_root()) {
  reference_composite <- file.path(
    project_root,
    "figures",
    "Figure_5",
    "output",
    "Figure_5_FinalFigure.png"
  )
  if (!file.exists(reference_composite)) {
    stop("The supplied Figure 5 reference composite is missing.", call. = FALSE)
  }
  if (!file.exists(result$file) || file.info(result$file)$size < 10000) {
    stop("Figure 5D PNG was not generated or is unexpectedly small.", call. = FALSE)
  }
  if (!identical(unname(read_png_dimensions(result$file)), c(825, 600))) {
    stop("Figure 5D PNG dimensions must match the 825 x 600 reference panel.", call. = FALSE)
  }
  invisible(reference_composite)
}

# - Figure generation
generate_figure_5_d <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_5", "output")) {
  project_root <- find_project_root()
  data <- prepare_figure_5_d_data(project_root)
  statistics <- calculate_figure_5_d_statistics(data)
  file <- file.path(output_dir, "Figure5D_generated.png")
  draw_figure_5_d(data, statistics, file)
  result <- list(
    file = file,
    data = data,
    statistics = statistics,
    source_files = attr(data, "source_files")
  )
  validate_figure_5_d_output(result, project_root)
  invisible(result)
}

if (sys.nframe() == 0) {
  result <- generate_figure_5_d()
  message("Generated Figure 5 panel D: ", normalizePath(result$file, winslash = "/"))
  print(result$statistics, row.names = FALSE)
}
