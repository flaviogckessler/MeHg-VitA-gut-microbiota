# Shared nonparametric taxon-abundance workflow for Figures S3-S5.
#
# The supplied final figures were assembled manually in Adobe Photoshop from
# individual 600 x 500 TIFF panels. The portable workflows generate only the
# corresponding individual PNG panels and never assemble the composites.

supplementary_abundance_required_packages <- c(
  "dplyr", "ggplot2", "ggpubr", "rstatix"
)
supplementary_abundance_missing_packages <-
  supplementary_abundance_required_packages[
    !vapply(
      supplementary_abundance_required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
if (length(supplementary_abundance_missing_packages) > 0L) {
  stop(
    "Install the packages required by Figures S3-S5: ",
    paste(supplementary_abundance_missing_packages, collapse = ", "),
    call. = FALSE
  )
}

supplementary_abundance_batch_ids <- c(
  batch1 = "batch01_iontorrent_dams-f0_f1-pnd21",
  batch2 = "batch02_illumina-miseq_f1-f2_pnd365"
)
supplementary_abundance_treatments <- c("CN", "VitA", "MeHg", "VitA+MeHg")

supplementary_abundance_panel_codes <- function(n) {
  encode <- function(index) {
    result <- character()
    while (index > 0L) {
      remainder <- (index - 1L) %% 26L
      result <- c(LETTERS[[remainder + 1L]], result)
      index <- (index - 1L) %/% 26L
    }
    paste0(result, collapse = "")
  }
  vapply(seq_len(n), encode, character(1))
}

find_supplementary_abundance_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata"),
      file.path(path, "results", "taxon_relative_abundance"),
      file.path(path, "supplementary_figures")
    )
    if (all(dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop(
        "Could not locate the repository root for supplementary abundance figures.",
        call. = FALSE
      )
    }
    path <- parent
  }
}

supplementary_abundance_source_files <- function(
    specification,
    project_root = find_supplementary_abundance_project_root()) {
  list(
    abundance = file.path(
      project_root,
      "results",
      "taxon_relative_abundance",
      paste0(specification$batch_id, "_taxon_relative_abundance.csv")
    ),
    metadata = file.path(
      project_root,
      "data",
      "metadata",
      specification$batch_id,
      "sample_manifest.csv"
    )
  )
}

supplementary_abundance_primary_flag <- function(values, filename) {
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

strip_supplementary_abundance_utf8_bom <- function(data) {
  if (length(names(data)) == 0L) return(data)

  first_header <- charToRaw(names(data)[[1L]])
  utf8_bom <- as.raw(c(239, 187, 191))
  if (length(first_header) >= 3L && identical(first_header[1:3], utf8_bom)) {
    names(data)[[1L]] <- rawToChar(first_header[-(1:3)])
  }
  data
}

validate_supplementary_abundance_specification <- function(specification) {
  required <- c(
    "figure_id", "batch_id", "generation", "panels", "treatment_levels",
    "treatment_colors", "expected_sample_counts", "expected_treatment_counts",
    "composite_filename", "composite_dimensions", "expected_abundance_md5",
    "expected_manifest_md5"
  )
  missing <- setdiff(required, names(specification))
  if (length(missing) > 0L) {
    stop(
      "Supplementary abundance specification is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!specification$batch_id %in% supplementary_abundance_batch_ids) {
    stop("Unexpected sequencing batch in figure specification.", call. = FALSE)
  }
  if (!identical(
    sort(specification$treatment_levels),
    sort(supplementary_abundance_treatments)
  )) {
    stop("Figure specification has invalid treatment levels.", call. = FALSE)
  }
  if (!identical(names(specification$treatment_colors), specification$treatment_levels)) {
    stop("Treatment colors must be named in plotting order.", call. = FALSE)
  }

  panels <- specification$panels
  panel_required <- c("panel", "sex", "label")
  if (!all(panel_required %in% names(panels)) || nrow(panels) == 0L) {
    stop("Figure specification has an invalid panel table.", call. = FALSE)
  }
  if (anyDuplicated(panels$panel) ||
      !identical(panels$panel, supplementary_abundance_panel_codes(nrow(panels)))) {
    stop("Panel IDs must be unique and sequential.", call. = FALSE)
  }
  if (anyDuplicated(paste(panels$sex, panels$label, sep = "|"))) {
    stop("Figure specification contains a duplicate sex-taxon panel.", call. = FALSE)
  }
  if (!all(panels$sex %in% names(specification$expected_sample_counts))) {
    stop("Panel table contains an unexpected sex.", call. = FALSE)
  }

  invisible(specification)
}

read_supplementary_abundance_sources <- function(
    specification,
    project_root = find_supplementary_abundance_project_root()) {
  validate_supplementary_abundance_specification(specification)
  files <- supplementary_abundance_source_files(specification, project_root)
  all_files <- unlist(files, use.names = FALSE)
  if (!all(file.exists(all_files))) {
    stop(
      "Missing source file(s): ",
      paste(all_files[!file.exists(all_files)], collapse = ", "),
      call. = FALSE
    )
  }

  abundance <- utils::read.csv(
    files$abundance,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  abundance <- strip_supplementary_abundance_utf8_bom(abundance)
  abundance_required <- c(
    "Num", "taxon_id", "tax_rank", "tax_name", "regex_match"
  )
  if (!all(abundance_required %in% names(abundance))) {
    stop("Taxon-abundance source has an invalid schema.", call. = FALSE)
  }

  metadata <- utils::read.csv(
    files$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata <- strip_supplementary_abundance_utf8_bom(metadata)
  metadata_required <- c(
    "sample_id", "batch_id", "generation", "treatment", "sex",
    "include_primary_analysis"
  )
  if (!all(metadata_required %in% names(metadata))) {
    stop("Sample manifest has an invalid schema.", call. = FALSE)
  }
  if (anyDuplicated(metadata$sample_id)) {
    stop("Sample manifest contains duplicate sample IDs.", call. = FALSE)
  }
  if (!all(metadata$batch_id == specification$batch_id)) {
    stop("Sample manifest contains an unexpected batch ID.", call. = FALSE)
  }
  if (!all(metadata$treatment %in% supplementary_abundance_treatments)) {
    stop("Sample manifest contains an unexpected treatment.", call. = FALSE)
  }
  if (!all(metadata$sex %in% c("Female", "Male"))) {
    stop("Sample manifest contains an unexpected sex.", call. = FALSE)
  }
  metadata$include_primary_analysis <- supplementary_abundance_primary_flag(
    metadata$include_primary_analysis,
    files$metadata
  )

  list(abundance = abundance, metadata = metadata, files = files)
}

resolve_supplementary_abundance_taxa <- function(abundance, labels) {
  normalized <- make.names(abundance$regex_match)
  resolved <- lapply(labels, function(label) {
    matches <- which(!is.na(normalized) & normalized == label)
    if (length(matches) != 1L) {
      stop(
        "Display label must resolve to exactly one canonical taxon row: ",
        label,
        call. = FALSE
      )
    }
    abundance$regex_match[[matches]]
  })
  unlist(resolved, use.names = FALSE)
}

prepare_supplementary_abundance_panel_data <- function(
    specification,
    sources = read_supplementary_abundance_sources(specification)) {
  panels <- specification$panels
  panels$taxon <- resolve_supplementary_abundance_taxa(
    sources$abundance,
    panels$label
  )
  metadata <- sources$metadata[
    sources$metadata$include_primary_analysis &
      sources$metadata$generation == specification$generation,
    ,
    drop = FALSE
  ]
  if (nrow(metadata) != sum(specification$expected_sample_counts)) {
    stop("Figure cohort has an unexpected total sample count.", call. = FALSE)
  }

  for (current_sex in names(specification$expected_sample_counts)) {
    sex_metadata <- metadata[metadata$sex == current_sex, , drop = FALSE]
    expected_n <- unname(specification$expected_sample_counts[[current_sex]])
    if (nrow(sex_metadata) != expected_n) {
      stop("Unexpected sample count for ", current_sex, ".", call. = FALSE)
    }
    observed_counts <- table(factor(
      sex_metadata$treatment,
      levels = specification$treatment_levels
    ))
    expected_counts <- specification$expected_treatment_counts[[current_sex]][
      specification$treatment_levels
    ]
    if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
      stop("Unexpected treatment counts for ", current_sex, ".", call. = FALSE)
    }
  }

  missing_samples <- setdiff(metadata$sample_id, names(sources$abundance))
  if (length(missing_samples) > 0L) {
    stop(
      "Manifest samples absent from the abundance source: ",
      paste(missing_samples, collapse = ", "),
      call. = FALSE
    )
  }

  result <- lapply(seq_len(nrow(panels)), function(index) {
    configuration <- panels[index, , drop = FALSE]
    panel_metadata <- metadata[
      metadata$sex == configuration$sex[[1L]],
      ,
      drop = FALSE
    ]
    taxon_row <- sources$abundance[
      sources$abundance$regex_match == configuration$taxon[[1L]],
      ,
      drop = FALSE
    ]
    if (nrow(taxon_row) != 1L) {
      stop("Resolved taxon row is not unique.", call. = FALSE)
    }
    values <- as.numeric(
      taxon_row[1L, panel_metadata$sample_id, drop = TRUE]
    )
    if (anyNA(values) || !all(is.finite(values)) || any(values < 0 | values > 1)) {
      stop(
        "Invalid relative abundance in panel ",
        configuration$panel[[1L]],
        ".",
        call. = FALSE
      )
    }

    data.frame(
      sample_id = panel_metadata$sample_id,
      treatment = factor(
        panel_metadata$treatment,
        levels = specification$treatment_levels
      ),
      relative_abundance = values,
      sex = configuration$sex[[1L]],
      taxon = configuration$taxon[[1L]],
      label = configuration$label[[1L]],
      panel = configuration$panel[[1L]],
      stringsAsFactors = FALSE
    )
  })
  names(result) <- panels$panel
  result
}

calculate_supplementary_abundance_statistics <- function(data) {
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

make_supplementary_abundance_panel <- function(data, specification) {
  statistics <- calculate_supplementary_abundance_statistics(data)
  plot <- ggpubr::ggboxplot(
    data,
    x = "treatment",
    y = "relative_abundance",
    fill = "treatment",
    outlier.shape = NA,
    palette = unname(specification$treatment_colors)
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
    ggplot2::labs(
      title = unique(data$label),
      subtitle = rstatix::get_test_label(
        statistics$omnibus,
        detailed = TRUE,
        type = "text"
      ),
      caption = rstatix::get_pwc_label(statistics$pairwise)
    )

  if (nrow(statistics$significant) > 0L) {
    plot <- plot + ggpubr::stat_pvalue_manual(
      statistics$significant,
      label = "p.adj.signif",
      tip.length = 0.01
    )
  }

  list(plot = plot, statistics = statistics)
}

generate_supplementary_abundance_figure <- function(
    specification,
    output_dir = file.path(
      find_supplementary_abundance_project_root(),
      "supplementary_figures",
      paste0("Figure_", specification$figure_id),
      "output"
    )) {
  validate_supplementary_abundance_specification(specification)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  panel_data <- prepare_supplementary_abundance_panel_data(specification)
  panels <- lapply(
    panel_data,
    make_supplementary_abundance_panel,
    specification = specification
  )
  plots <- lapply(panels, `[[`, "plot")
  statistics <- lapply(panels, `[[`, "statistics")
  files <- setNames(
    file.path(
      output_dir,
      paste0(
        "Figure",
        specification$figure_id,
        specification$panels$panel,
        "_generated.png"
      )
    ),
    specification$panels$panel
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

check_supplementary_abundance_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

check_supplementary_abundance_numeric <- function(
    observed,
    expected,
    tolerance,
    message) {
  check_supplementary_abundance_true(
    length(observed) == length(expected) &&
      all(is.finite(observed)) &&
      isTRUE(all.equal(
        as.numeric(observed),
        as.numeric(expected),
        tolerance = tolerance,
        check.attributes = FALSE
      )),
    message
  )
}

read_supplementary_abundance_png_dimensions <- function(filename) {
  header <- readBin(filename, what = "raw", n = 24L)
  png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(header) != 24L || !identical(header[1:8], png_signature)) {
    stop("Invalid PNG file: ", filename, call. = FALSE)
  }
  uint32 <- function(bytes) sum(as.integer(bytes) * 256^(3:0))
  c(
    width = uint32(header[17:20]),
    height = uint32(header[21:24])
  )
}

validate_supplementary_abundance_figure <- function(
    specification,
    generator) {
  validate_supplementary_abundance_specification(specification)
  project_root <- find_supplementary_abundance_project_root()
  sources <- read_supplementary_abundance_sources(specification, project_root)
  observed_md5 <- unname(tools::md5sum(unlist(sources$files, use.names = FALSE)))
  expected_md5 <- tolower(c(
    specification$expected_abundance_md5,
    specification$expected_manifest_md5
  ))
  check_supplementary_abundance_true(
    identical(tolower(observed_md5), expected_md5),
    paste0("Figure ", specification$figure_id, " source checksum changed.")
  )

  panel_data <- prepare_supplementary_abundance_panel_data(
    specification,
    sources
  )
  check_supplementary_abundance_true(
    identical(names(panel_data), specification$panels$panel),
    "Prepared panel order changed."
  )
  check_supplementary_abundance_true(
    all(vapply(panel_data, function(data) {
      nrow(data) == specification$expected_sample_counts[[unique(data$sex)]] &&
        length(unique(data$sample_id)) == nrow(data) &&
        all(is.finite(data$relative_abundance)) &&
        all(data$relative_abundance >= 0 & data$relative_abundance <= 1) &&
        identical(make.names(unique(data$taxon)), unique(data$label))
    }, logical(1))),
    "Panel data, sample counts, or taxon labels changed."
  )

  statistics <- lapply(panel_data, calculate_supplementary_abundance_statistics)
  check_supplementary_abundance_true(
    all(vapply(statistics, function(result) {
      identical(result$omnibus$method, "Kruskal-Wallis") &&
        result$omnibus$df[[1L]] == 3L &&
        result$omnibus$p[[1L]] <= 0.05 &&
        nrow(result$pairwise) == 6L &&
        all(is.finite(result$pairwise$p)) &&
        all(is.finite(result$pairwise$p.adj)) &&
        isTRUE(all.equal(
          result$pairwise$p.adj,
          stats::p.adjust(result$pairwise$p, method = "BH"),
          tolerance = 1e-14,
          check.attributes = FALSE
        ))
    }, logical(1))),
    "Nonparametric test coverage or FDR adjustment changed."
  )

  independent_omnibus <- vapply(panel_data, function(data) {
    unname(stats::kruskal.test(
      relative_abundance ~ treatment,
      data = data
    )$statistic)
  }, numeric(1))
  observed_omnibus <- vapply(statistics, function(result) {
    result$omnibus$statistic[[1L]]
  }, numeric(1))
  check_supplementary_abundance_numeric(
    observed_omnibus,
    independent_omnibus,
    tolerance = 1e-13,
    "rstatix and base-R Kruskal-Wallis statistics disagree."
  )

  if (!is.null(specification$expected_omnibus_statistics)) {
    check_supplementary_abundance_numeric(
      observed_omnibus,
      specification$expected_omnibus_statistics,
      tolerance = 1e-12,
      "Expected Kruskal-Wallis statistics changed."
    )
  }
  if (!is.null(specification$expected_omnibus_p)) {
    observed_p <- vapply(statistics, function(result) {
      result$omnibus$p[[1L]]
    }, numeric(1))
    check_supplementary_abundance_numeric(
      observed_p,
      specification$expected_omnibus_p,
      tolerance = 1e-12,
      "Expected Kruskal-Wallis p-values changed."
    )
  }
  if (!is.null(specification$expected_significant_counts)) {
    observed_counts <- vapply(statistics, function(result) {
      nrow(result$significant)
    }, integer(1))
    check_supplementary_abundance_true(
      identical(
        as.integer(observed_counts),
        as.integer(specification$expected_significant_counts)
      ),
      "Expected significant Dunn-comparison counts changed."
    )
  }
  if (!is.null(specification$expected_significant_signatures)) {
    observed_signatures <- vapply(statistics, function(result) {
      significant <- result$significant
      if (nrow(significant) == 0L) return("")
      paste(
        paste(significant$group1, significant$group2, sep = "|"),
        collapse = ","
      )
    }, character(1))
    check_supplementary_abundance_true(
      identical(
        unname(observed_signatures),
        unname(specification$expected_significant_signatures)
      ),
      "Expected significant Dunn-comparison pairs changed."
    )
  }

  validation_dir <- file.path(
    tempdir(),
    paste0("figure_", tolower(specification$figure_id), "_validation")
  )
  rendered <- generator(validation_dir)
  expected_names <- paste0(
    "Figure",
    specification$figure_id,
    specification$panels$panel,
    "_generated.png"
  )
  check_supplementary_abundance_true(
    identical(basename(unname(rendered$files)), expected_names) &&
      all(file.exists(rendered$files)),
    "Generated panel files are missing or misnamed."
  )
  check_supplementary_abundance_true(
    all(vapply(rendered$files, function(filename) {
      identical(
        unname(read_supplementary_abundance_png_dimensions(filename)),
        c(600, 500)
      )
    }, logical(1))),
    "Generated panel PNGs must be 600 x 500 pixels."
  )
  check_supplementary_abundance_true(
    all(file.info(rendered$files)$size > 5000 &
        file.info(rendered$files)$size < 100000) &&
      all(vapply(rendered$plots, inherits, logical(1), what = "ggplot")),
    "Generated panel PNGs or ggplot objects are invalid."
  )

  committed_files <- file.path(
    project_root,
    "supplementary_figures",
    paste0("Figure_", specification$figure_id),
    "output",
    expected_names
  )
  check_supplementary_abundance_true(
    all(file.exists(committed_files)) &&
      all(vapply(committed_files, function(filename) {
        identical(
          unname(read_supplementary_abundance_png_dimensions(filename)),
          c(600, 500)
        )
      }, logical(1))),
    "Committed generated panels are missing or have invalid dimensions."
  )

  composite <- file.path(
    project_root,
    "supplementary_figures",
    paste0("Figure_", specification$figure_id),
    "output",
    specification$composite_filename
  )
  check_supplementary_abundance_true(
    file.exists(composite) &&
      identical(
        unname(read_supplementary_abundance_png_dimensions(composite)),
        as.numeric(specification$composite_dimensions)
      ),
    "Supplied reference composite is missing or has changed dimensions."
  )

  invisible(list(
    panel_data = panel_data,
    statistics = statistics,
    rendered = rendered
  ))
}
