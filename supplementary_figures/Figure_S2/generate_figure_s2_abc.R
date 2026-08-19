# Reproducible Figure S2 panels A-C.
#
# The supplied Figure_S2_FinalFigure.png was assembled manually in Adobe
# Photoshop. This script generates only the three individual Dams (F0)
# relative-abundance boxplots used in that layout. The legacy workflow is
# preserved in scripts/legacy_scripts/Abundancia_boxplots.R and is not sourced
# at runtime.

figure_s2_required_packages <- c("dplyr", "ggplot2", "ggpubr", "rstatix")
figure_s2_missing_packages <- figure_s2_required_packages[
  !vapply(
    figure_s2_required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(figure_s2_missing_packages) > 0L) {
  stop(
    "Install the packages required by Figure S2: ",
    paste(figure_s2_missing_packages, collapse = ", "),
    call. = FALSE
  )
}

figure_s2_batch_id <- "batch01_iontorrent_dams-f0_f1-pnd21"
figure_s2_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
figure_s2_colors <- c(
  CN = "#98F5FF",
  VitA = "#B3EE3A",
  MeHg = "#CD2626",
  `VitA+MeHg` = "#4D4D4D"
)
figure_s2_panel_configurations <- list(
  A = list(
    taxon = "g__[Eubacterium] ventriosum group",
    label = "g__.Eubacterium..ventriosum.group"
  ),
  B = list(
    taxon = "g__[Ruminococcus] gauvreauii group",
    label = "g__.Ruminococcus..gauvreauii.group"
  ),
  C = list(
    taxon = "g__Oribacterium",
    label = "g__Oribacterium"
  )
)

find_figure_s2_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(
        path,
        "results",
        "taxon_relative_abundance",
        paste0(figure_s2_batch_id, "_taxon_relative_abundance.csv")
      ),
      file.path(
        path,
        "data",
        "metadata",
        figure_s2_batch_id,
        "sample_manifest.csv"
      ),
      file.path(path, "supplementary_figures", "Figure_S2")
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the repository root for Figure S2.", call. = FALSE)
    }
    path <- parent
  }
}

figure_s2_source_files <- function(
    project_root = find_figure_s2_project_root()) {
  list(
    abundance = file.path(
      project_root,
      "results",
      "taxon_relative_abundance",
      paste0(figure_s2_batch_id, "_taxon_relative_abundance.csv")
    ),
    metadata = file.path(
      project_root,
      "data",
      "metadata",
      figure_s2_batch_id,
      "sample_manifest.csv"
    )
  )
}

figure_s2_primary_flag <- function(values, filename) {
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

strip_figure_s2_utf8_bom <- function(data) {
  if (length(names(data)) == 0L) return(data)

  first_header <- charToRaw(names(data)[[1L]])
  utf8_bom <- as.raw(c(239, 187, 191))
  if (length(first_header) >= 3L && identical(first_header[1:3], utf8_bom)) {
    names(data)[[1L]] <- rawToChar(first_header[-(1:3)])
  }
  data
}

read_figure_s2_sources <- function(
    project_root = find_figure_s2_project_root()) {
  files <- figure_s2_source_files(project_root)
  missing_files <- unlist(files, use.names = FALSE)[
    !file.exists(unlist(files, use.names = FALSE))
  ]
  if (length(missing_files) > 0L) {
    stop(
      "Missing Figure S2 source file(s): ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }

  abundance <- utils::read.csv(
    files$abundance,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  abundance <- strip_figure_s2_utf8_bom(abundance)
  abundance_required <- c(
    "Num", "taxon_id", "tax_rank", "tax_name", "regex_match"
  )
  if (!all(abundance_required %in% names(abundance))) {
    stop("The Figure S2 abundance source has an invalid schema.", call. = FALSE)
  }
  metadata <- utils::read.csv(
    files$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata <- strip_figure_s2_utf8_bom(metadata)
  metadata_required <- c(
    "sample_id", "batch_id", "generation", "treatment", "sex",
    "include_primary_analysis"
  )
  if (!all(metadata_required %in% names(metadata))) {
    stop("The Figure S2 sample manifest has an invalid schema.", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("The Figure S2 sample manifest contains duplicate IDs.", call. = FALSE)
  }
  if (!all(metadata$batch_id == figure_s2_batch_id)) {
    stop("The Figure S2 manifest contains an unexpected batch ID.", call. = FALSE)
  }
  if (!all(metadata$treatment %in% figure_s2_treatment_levels)) {
    stop("The Figure S2 manifest contains an unexpected treatment.", call. = FALSE)
  }
  metadata$include_primary_analysis <- figure_s2_primary_flag(
    metadata$include_primary_analysis,
    files$metadata
  )

  list(abundance = abundance, metadata = metadata, files = files)
}

prepare_figure_s2_panel_data <- function(
    sources = read_figure_s2_sources()) {
  configurations <- figure_s2_panel_configurations
  required_taxa <- vapply(configurations, `[[`, character(1), "taxon")
  missing_taxa <- setdiff(required_taxa, sources$abundance$regex_match)
  if (length(missing_taxa) > 0L) {
    stop(
      "Figure S2 taxa absent from the abundance source: ",
      paste(missing_taxa, collapse = ", "),
      call. = FALSE
    )
  }
  selected_taxon_counts <- table(
    factor(sources$abundance$regex_match, levels = required_taxa),
    useNA = "no"
  )
  if (any(selected_taxon_counts != 1L)) {
    stop("Each Figure S2 taxon must occur exactly once.", call. = FALSE)
  }

  metadata <- sources$metadata[
    sources$metadata$include_primary_analysis &
      sources$metadata$generation == "F0",
    ,
    drop = FALSE
  ]
  if (nrow(metadata) != 18L || !all(metadata$sex == "Female")) {
    stop("Figure S2 must contain the 18 primary-analysis F0 dams.", call. = FALSE)
  }
  missing_samples <- setdiff(metadata$sample_id, names(sources$abundance))
  if (length(missing_samples) > 0L) {
    stop(
      "Figure S2 samples absent from the abundance source: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  panels <- lapply(configurations, function(configuration) {
    taxon_row <- sources$abundance[
      sources$abundance$regex_match == configuration$taxon,
      ,
      drop = FALSE
    ]
    values <- as.numeric(taxon_row[1L, metadata$sample_id, drop = TRUE])
    if (anyNA(values) || !all(is.finite(values)) || any(values < 0 | values > 1)) {
      stop(
        "Figure S2 contains an invalid relative abundance for ",
        configuration$taxon,
        ".",
        call. = FALSE
      )
    }

    data.frame(
      sample_id = metadata$sample_id,
      treatment = factor(
        metadata$treatment,
        levels = figure_s2_treatment_levels
      ),
      relative_abundance = values,
      taxon = configuration$taxon,
      stringsAsFactors = FALSE
    )
  })
  names(panels) <- names(configurations)
  panels
}

calculate_figure_s2_statistics <- function(data) {
  omnibus <- rstatix::kruskal_test(
    data,
    relative_abundance ~ treatment
  )
  pairwise <- rstatix::dunn_test(
    data,
    relative_abundance ~ treatment,
    p.adjust.method = "fdr"
  )
  positioned <- rstatix::add_y_position(
    pairwise,
    formula = relative_abundance ~ treatment,
    data = data
  )
  significant <- positioned[
    !is.na(positioned$p.adj) & positioned$p.adj <= 0.05,
    ,
    drop = FALSE
  ]

  list(
    omnibus = omnibus,
    pairwise = pairwise,
    significant = significant
  )
}

make_figure_s2_panel <- function(data, configuration) {
  statistics <- calculate_figure_s2_statistics(data)

  plot <- ggpubr::ggboxplot(
    data,
    x = "treatment",
    y = "relative_abundance",
    fill = "treatment",
    outlier.shape = NA,
    palette = unname(figure_s2_colors)
  ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.direction = "vertical",
      legend.position = "right",
      panel.background = ggplot2::element_rect(
        fill = "gray93",
        colour = "grey",
        linetype = "solid"
      ),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.5,
        linetype = "solid",
        colour = "white"
      ),
      plot.margin = ggplot2::margin(
        t = 4.06,
        r = 5.5,
        b = 6.94,
        l = 5.5,
        unit = "pt"
      )
    ) +
    ggplot2::xlab("Groups") +
    ggplot2::ylab("Relative Abundance") +
    ggpubr::stat_pvalue_manual(
      statistics$significant,
      label = "p.adj.signif",
      tip.length = 0.01
    ) +
    ggplot2::labs(
      title = configuration$label,
      subtitle = rstatix::get_test_label(
        statistics$omnibus,
        detailed = TRUE,
        type = "text"
      ),
      caption = rstatix::get_pwc_label(statistics$significant)
    )

  list(plot = plot, statistics = statistics)
}

generate_figure_s2_abc <- function(
    output_dir = file.path(
      find_figure_s2_project_root(),
      "supplementary_figures",
      "Figure_S2",
      "output"
    )) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  panel_data <- prepare_figure_s2_panel_data()
  panels <- Map(
    make_figure_s2_panel,
    panel_data,
    figure_s2_panel_configurations
  )
  plots <- lapply(panels, `[[`, "plot")
  statistics <- lapply(panels, `[[`, "statistics")
  files <- setNames(
    file.path(
      output_dir,
      paste0("FigureS2", names(figure_s2_panel_configurations), "_generated.png")
    ),
    names(figure_s2_panel_configurations)
  )

  Map(function(filename, plot) {
    ggplot2::ggsave(
      filename = filename,
      plot = plot,
      device = "png",
      width = 6,
      height = 5,
      units = "in",
      dpi = 100
    )
  }, files, plots)

  invisible(list(
    files = files,
    plots = plots,
    data = panel_data,
    statistics = statistics
  ))
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s2_abc()
  message(
    "Generated Figure S2 panels A-C in: ",
    normalizePath(dirname(result$files[[1L]]), winslash = "/")
  )
}
