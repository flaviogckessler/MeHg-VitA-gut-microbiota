# Reproducible Figure 8 panel E: KO z-score heatmap for KEGG map00480.
#
# This is a portable base-R refactor of the heatmap.2() workflow near line 981
# of scripts/legacy_scripts/Picrust_mestrado_Alberto.R. The legacy script is
# retained for provenance and is not a runtime dependency. Panel E reads the
# clean F1-PND21/map00480 table generated from canonical PICRUSt2 KO abundance
# by scripts/04_functional_prediction/generate_kegg_pathway_ko_z_scores.R.

# Load Packages ----
suppressPackageStartupMessages({
  library(graphics)
  library(grDevices)
  library(stats)
  library(utils)
})

figure_8_e_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")

# The legacy heatmap used hierarchical clustering without column clustering.
# All profiles already have row mean zero, so dendrogram rotations can change
# across platforms when branches are reordered using numerically tied means.
# This checkpoint records the supplied reference panel's top-to-bottom legal
# rotation of the same complete-linkage dendrogram.
figure_8_e_reference_row_order <- c(
  "K00681", "K01919", "K01920", "K00799", "K00032", "K07751",
  "K01581", "K01256", "K00383", "K00033", "K00036", "K00797",
  "K01460", "K00432", "K01255", "K01469", "K00031", "K01270"
)

find_figure_8_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(
        path, "results", "functional_prediction",
        "batch01_iontorrent_dams-f0_f1-pnd21",
        "kegg_pathway_ko_z_scores",
        "f1_pnd21_ko00480_ko_z_scores.csv"
      ),
      file.path(path, "figures", "Figure_8"),
      file.path(
        path, "scripts", "04_functional_prediction",
        "generate_kegg_pathway_ko_z_scores.R"
      )
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the New_repo project root.", call. = FALSE)
    }
    path <- parent
  }
}

read_figure_8_e_data <- function(project_root = find_figure_8_project_root()) {
  source_file <- file.path(
    project_root,
    "results", "functional_prediction",
    "batch01_iontorrent_dams-f0_f1-pnd21",
    "kegg_pathway_ko_z_scores",
    "f1_pnd21_ko00480_ko_z_scores.csv"
  )
  data <- utils::read.csv(
    source_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )

  expected_columns <- c("ko_id", figure_8_e_treatment_levels)
  if (!identical(names(data), expected_columns)) {
    stop(
      "The Figure 8E source columns must be: ",
      paste(expected_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(data$ko_id) || any(!grepl("^K[0-9]{5}$", data$ko_id))) {
    stop("The Figure 8E source contains invalid or duplicated KO IDs.", call. = FALSE)
  }
  if (!setequal(data$ko_id, figure_8_e_reference_row_order)) {
    stop(
      "The measured map00480 KO set differs from the supplied Figure 8E reference.",
      call. = FALSE
    )
  }

  z_scores <- as.matrix(data[, figure_8_e_treatment_levels, drop = FALSE])
  storage.mode(z_scores) <- "double"
  rownames(z_scores) <- data$ko_id
  if (any(!is.finite(z_scores))) {
    stop("The Figure 8E source contains non-finite z-scores.", call. = FALSE)
  }
  if (
    any(abs(rowMeans(z_scores)) > 1e-12) ||
      any(abs(apply(z_scores, 1L, stats::sd) - 1) > 1e-12)
  ) {
    stop("Figure 8E rows must have mean zero and unit sample SD.", call. = FALSE)
  }

  attr(z_scores, "source_file") <- source_file
  z_scores
}

orient_figure_8_e_dendrogram <- function(dendrogram, row_order) {
  target_position <- setNames(seq_along(row_order), row_order)

  orient_node <- function(node) {
    if (is.leaf(node)) return(node)
    children <- lapply(node, orient_node)
    first_position <- vapply(children, function(child) {
      min(unname(target_position[labels(child)]))
    }, numeric(1))
    node[] <- children[order(first_position)]
    node
  }

  oriented <- orient_node(dendrogram)
  if (!identical(labels(oriented), row_order)) {
    stop(
      "The reference KO order is not a legal rotation of the source-data dendrogram.",
      call. = FALSE
    )
  }
  oriented
}

draw_figure_8_e_dendrogram <- function(
    dendrogram,
    row_centers,
    leaf_x = 14,
    root_x = 0) {
  label_y <- setNames(row_centers, labels(dendrogram))
  maximum_height <- attr(dendrogram, "height")

  draw_node <- function(node) {
    if (is.leaf(node)) {
      return(c(x = leaf_x, y = unname(label_y[attr(node, "label")])))
    }

    child_coordinates <- lapply(node, draw_node)
    node_x <- leaf_x -
      attr(node, "height") / maximum_height * (leaf_x - root_x)
    child_y <- vapply(child_coordinates, `[[`, numeric(1), "y")
    graphics::segments(
      node_x, min(child_y), node_x, max(child_y),
      col = "black", lwd = 0.80
    )
    for (child in child_coordinates) {
      graphics::segments(
        node_x, child[["y"]], child[["x"]], child[["y"]],
        col = "black", lwd = 0.80
      )
    }
    c(x = node_x, y = mean(range(child_y)))
  }

  invisible(draw_node(dendrogram))
}

draw_figure_8_e <- function(z_scores, filename) {
  row_order <- figure_8_e_reference_row_order
  z_scores <- z_scores[row_order, figure_8_e_treatment_levels, drop = FALSE]
  row_count <- nrow(z_scores)

  dendrogram <- as.dendrogram(stats::hclust(stats::dist(z_scores)))
  dendrogram <- orient_figure_8_e_dendrogram(dendrogram, row_order)

  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  png_type <- if (capabilities("cairo")) "cairo" else "windows"
  grDevices::png(
    filename,
    width = 600,
    height = 1800,
    units = "px",
    res = 150,
    pointsize = 12,
    bg = "white",
    type = png_type
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  graphics::par(
    mai = c(0, 0, 0, 0),
    family = "sans",
    xaxs = "i",
    yaxs = "i"
  )
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 600),
    ylim = c(1800, 0),
    xaxs = "i",
    yaxs = "i"
  )

  # Pixel coordinates reconstructed from the supplied 600 x 1800 reference.
  heatmap_left <- 14
  column_edges <- c(14, 111, 207, 304, 400)
  row_edges <- c(
    43, 130, 216, 303, 389, 476, 562, 649, 735, 822,
    908, 995, 1081, 1168, 1254, 1341, 1427, 1514, 1600
  )
  row_centers <- (row_edges[-length(row_edges)] + row_edges[-1L]) / 2

  palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(15L)
  color_limit <- max(abs(z_scores))
  breaks <- seq(-color_limit, color_limit, length.out = length(palette) + 1L)
  color_index <- findInterval(
    z_scores,
    breaks,
    all.inside = TRUE,
    left.open = FALSE
  )
  dim(color_index) <- dim(z_scores)

  for (row_index in seq_len(row_count)) {
    for (column_index in seq_along(figure_8_e_treatment_levels)) {
      graphics::rect(
        xleft = column_edges[column_index],
        ybottom = row_edges[row_index],
        xright = column_edges[column_index + 1L],
        ytop = row_edges[row_index + 1L],
        col = palette[color_index[row_index, column_index]],
        border = NA
      )
    }
  }

  draw_figure_8_e_dendrogram(
    dendrogram = dendrogram,
    row_centers = row_centers,
    leaf_x = heatmap_left,
    root_x = 0
  )
  graphics::text(
    x = 425,
    y = row_centers + 0.5,
    labels = row_order,
    adj = c(0, 0.5),
    cex = 0.83,
    col = "black",
    xpd = NA
  )
  graphics::text(
    x = c(55, 152, 248, 345),
    y = c(1626, 1627, 1628, 1630),
    labels = figure_8_e_treatment_levels,
    srt = 90,
    adj = c(1, 0.5),
    cex = 0.83,
    col = "black",
    xpd = NA
  )

  invisible(filename)
}

read_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)
  header <- readBin(connection, what = "raw", n = 24L)
  signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(header) < 24L || !identical(header[1:8], signature)) {
    stop("Invalid PNG signature: ", filename, call. = FALSE)
  }
  unsigned_32 <- function(bytes) {
    sum(as.numeric(as.integer(bytes)) * 256^(3:0))
  }
  c(
    width = unsigned_32(header[17:20]),
    height = unsigned_32(header[21:24])
  )
}

read_tiff_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)

  unsigned_integer <- function(bytes, endian) {
    values <- as.numeric(as.integer(bytes))
    powers <- if (identical(endian, "little")) {
      seq_along(values) - 1L
    } else {
      rev(seq_along(values) - 1L)
    }
    sum(values * 256^powers)
  }

  header <- readBin(connection, what = "raw", n = 8L)
  if (length(header) < 8L) {
    stop("Invalid TIFF header: ", filename, call. = FALSE)
  }
  byte_order <- rawToChar(header[1:2])
  endian <- switch(
    byte_order,
    II = "little",
    MM = "big",
    stop("Invalid TIFF byte order: ", filename, call. = FALSE)
  )
  if (unsigned_integer(header[3:4], endian) != 42) {
    stop("Invalid TIFF signature: ", filename, call. = FALSE)
  }

  ifd_offset <- unsigned_integer(header[5:8], endian)
  seek(connection, where = ifd_offset, origin = "start")
  entry_count_bytes <- readBin(connection, what = "raw", n = 2L)
  if (length(entry_count_bytes) != 2L) {
    stop("Invalid TIFF image-file directory: ", filename, call. = FALSE)
  }
  entry_count <- unsigned_integer(entry_count_bytes, endian)
  dimensions <- c(width = NA_real_, height = NA_real_)

  for (entry_index in seq_len(entry_count)) {
    entry <- readBin(connection, what = "raw", n = 12L)
    if (length(entry) != 12L) {
      stop("Truncated TIFF image-file directory: ", filename, call. = FALSE)
    }
    tag <- unsigned_integer(entry[1:2], endian)
    field_type <- unsigned_integer(entry[3:4], endian)
    value_count <- unsigned_integer(entry[5:8], endian)
    if (tag %in% c(256, 257) && value_count == 1 && field_type %in% c(3, 4)) {
      value_bytes <- if (field_type == 3) entry[9:10] else entry[9:12]
      dimension_name <- if (tag == 256) "width" else "height"
      dimensions[dimension_name] <- unsigned_integer(value_bytes, endian)
    }
  }

  if (any(!is.finite(dimensions))) {
    stop("TIFF width or height tag is missing: ", filename, call. = FALSE)
  }
  dimensions
}

validate_figure_8_e_output <- function(result, project_root) {
  reference_panel <- file.path(
    project_root, "figures", "Figure_8", "output", "Figure_7_heatmap.tiff"
  )
  reference_composite <- file.path(
    project_root, "figures", "Figure_8", "output", "Figure_8_FinalFigure.png"
  )
  if (!file.exists(result$file) || file.info(result$file)$size < 10000) {
    stop("Figure 8E PNG was not generated or is unexpectedly small.", call. = FALSE)
  }
  if (!identical(unname(read_png_dimensions(result$file)), c(600, 1800))) {
    stop("Figure 8E PNG dimensions must be 600 x 1800.", call. = FALSE)
  }

  missing_references <- c(reference_panel, reference_composite)[
    !file.exists(c(reference_panel, reference_composite))
  ]
  if (length(missing_references) > 0L) {
    message(
      "Generated Figure 8E successfully. Optional visual-reference validation ",
      "was skipped because the following file(s) are unavailable: ",
      paste(basename(missing_references), collapse = ", "),
      "."
    )
    return(invisible(TRUE))
  }

  if (!identical(unname(read_tiff_dimensions(reference_panel)), c(600, 1800))) {
    stop("The supplied standalone Figure 8E reference dimensions changed.", call. = FALSE)
  }
  if (!identical(unname(read_png_dimensions(reference_composite)), c(3449, 2915))) {
    stop("The supplied final Figure 8 reference dimensions changed.", call. = FALSE)
  }
  invisible(TRUE)
}

generate_figure_8_e <- function(
    output_dir = file.path(
      find_figure_8_project_root(), "figures", "Figure_8", "output"
    )) {
  project_root <- find_figure_8_project_root()
  z_scores <- read_figure_8_e_data(project_root)
  output_file <- file.path(output_dir, "Figure8E_generated.png")
  draw_figure_8_e(z_scores, output_file)
  result <- list(
    file = output_file,
    data = z_scores,
    row_order = figure_8_e_reference_row_order,
    source_file = attr(z_scores, "source_file")
  )
  validate_figure_8_e_output(result, project_root)
  invisible(result)
}

if (sys.nframe() == 0L) {
  result <- generate_figure_8_e()
  message(
    "Generated Figure 8 panel E: ",
    normalizePath(result$file, winslash = "/")
  )
}
