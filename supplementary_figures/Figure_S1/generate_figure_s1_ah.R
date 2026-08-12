# Reproducible Figure S1 panels A-H.
#
# Figure S1 complements Figure 3 panels D-G. Panels A-D show Observed richness
# and panels E-H show InvSimpson diversity for the same four cohorts. The
# supplied Figure_S1_FinalFigure.png was assembled manually in Adobe Photoshop;
# this script generates the eight individual PNG panels used for that layout.

figure_s1_required_packages <- c("dplyr", "ggplot2", "ggpubr", "rstatix")
figure_s1_missing_packages <- figure_s1_required_packages[
  !vapply(
    figure_s1_required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(figure_s1_missing_packages) > 0L) {
  stop(
    "Install the packages required by Figure S1: ",
    paste(figure_s1_missing_packages, collapse = ", "),
    call. = FALSE
  )
}

figure_s1_batch_ids <- c(
  batch1 = "batch01_iontorrent_dams-f0_f1-pnd21",
  batch2 = "batch02_illumina-miseq_f1-f2_pnd365"
)
figure_s1_treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
figure_s1_sex_levels <- c("Female", "Male")
figure_s1_colors <- c(
  CN = "#98F5FF",
  VitA = "#B3EE3A",
  MeHg = "#CD2626",
  `VitA+MeHg` = "#4D4D4D"
)

figure_s1_cohort_configurations <- list(
  DAMS_F0 = list(
    label = "Dams (F0)",
    batch_id = unname(figure_s1_batch_ids[["batch1"]]),
    generation = "F0",
    expected_sample_count = 18L,
    expected_sexes = "Female",
    width = 4
  ),
  F1_PND21 = list(
    label = "F1-PND21",
    batch_id = unname(figure_s1_batch_ids[["batch1"]]),
    generation = "F1",
    expected_sample_count = 34L,
    expected_sexes = figure_s1_sex_levels,
    width = 8
  ),
  F1_PND365 = list(
    label = "F1-PND365",
    batch_id = unname(figure_s1_batch_ids[["batch2"]]),
    generation = "F1",
    expected_sample_count = 24L,
    expected_sexes = figure_s1_sex_levels,
    width = 8
  ),
  F2_PND365 = list(
    label = "F2-PND365",
    batch_id = unname(figure_s1_batch_ids[["batch2"]]),
    generation = "F2",
    expected_sample_count = 24L,
    expected_sexes = figure_s1_sex_levels,
    width = 8
  )
)

figure_s1_panel_configurations <- list(
  A = list(metric = "Observed", cohort = "DAMS_F0", method = "parametric"),
  B = list(metric = "Observed", cohort = "F1_PND21", method = "nonparametric"),
  C = list(metric = "Observed", cohort = "F1_PND365", method = "nonparametric"),
  D = list(metric = "Observed", cohort = "F2_PND365", method = "parametric"),
  E = list(metric = "InvSimpson", cohort = "DAMS_F0", method = "parametric"),
  F = list(metric = "InvSimpson", cohort = "F1_PND21", method = "parametric"),
  G = list(metric = "InvSimpson", cohort = "F1_PND365", method = "parametric"),
  H = list(metric = "InvSimpson", cohort = "F2_PND365", method = "parametric")
)

find_figure_s1_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata", "sample_manifest_all_batches.csv"),
      file.path(path, "results", "alpha_diversity"),
      file.path(path, "supplementary_figures", "Figure_S1")
    )
    if (all(file.exists(markers) | dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the repository root for Figure S1.", call. = FALSE)
    }
    path <- parent
  }
}

figure_s1_source_files <- function(project_root = find_figure_s1_project_root()) {
  list(
    alpha = setNames(
      file.path(
        project_root,
        "results",
        "alpha_diversity",
        c(
          paste0(figure_s1_batch_ids[["batch1"]], "_alpha_diversity_with_metadata.csv"),
          paste0(figure_s1_batch_ids[["batch2"]], "_alpha_diversity_with_metadata.csv")
        )
      ),
      unname(figure_s1_batch_ids)
    ),
    metadata = setNames(
      file.path(
        project_root,
        "data",
        "metadata",
        unname(figure_s1_batch_ids),
        "sample_manifest.csv"
      ),
      unname(figure_s1_batch_ids)
    )
  )
}

figure_s1_primary_flag <- function(values, filename) {
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

strip_figure_s1_utf8_bom <- function(data) {
  if (length(names(data)) == 0L) return(data)

  first_header <- charToRaw(names(data)[[1]])
  utf8_bom <- as.raw(c(239, 187, 191))
  if (length(first_header) >= 3L && identical(first_header[1:3], utf8_bom)) {
    names(data)[[1]] <- rawToChar(first_header[-(1:3)])
  }
  data
}

read_figure_s1_alpha_source <- function(filename, sample_column) {
  data <- utils::read.csv(
    filename,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  data <- strip_figure_s1_utf8_bom(data)
  required <- c(sample_column, "Observed", "InvSimpson")
  if (!all(required %in% names(data))) {
    stop("Alpha-diversity source is missing required columns: ", filename, call. = FALSE)
  }

  result <- data.frame(
    sample_id = trimws(as.character(data[[sample_column]])),
    Observed = as.numeric(data[["Observed"]]),
    InvSimpson = as.numeric(data[["InvSimpson"]]),
    stringsAsFactors = FALSE
  )
  if (anyNA(result$sample_id) || any(result$sample_id == "")) {
    stop("Alpha-diversity source contains an invalid sample ID: ", filename, call. = FALSE)
  }
  if (anyDuplicated(result$sample_id)) {
    stop("Alpha-diversity source contains duplicated sample IDs: ", filename, call. = FALSE)
  }
  if (!all(is.finite(result$Observed)) || !all(is.finite(result$InvSimpson))) {
    stop("Alpha-diversity source contains non-finite values: ", filename, call. = FALSE)
  }
  if (any(result$Observed <= 0) || any(result$InvSimpson < 1)) {
    stop("Alpha-diversity values fall outside their valid range: ", filename, call. = FALSE)
  }

  result
}

read_figure_s1_manifest <- function(filename, expected_batch_id) {
  data <- utils::read.csv(
    filename,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  data <- strip_figure_s1_utf8_bom(data)
  required <- c(
    "sample_id", "batch_id", "generation", "treatment", "sex",
    "include_primary_analysis"
  )
  if (!all(required %in% names(data))) {
    stop("Sample manifest is missing required columns: ", filename, call. = FALSE)
  }
  if (anyDuplicated(data$sample_id)) {
    stop("Sample manifest contains duplicated sample IDs: ", filename, call. = FALSE)
  }
  if (!all(data$batch_id == expected_batch_id)) {
    stop("Sample manifest contains an unexpected batch ID: ", filename, call. = FALSE)
  }
  if (!all(data$treatment %in% figure_s1_treatment_levels)) {
    stop("Sample manifest contains an unexpected treatment: ", filename, call. = FALSE)
  }
  if (!all(data$sex %in% figure_s1_sex_levels)) {
    stop("Sample manifest contains an unexpected sex: ", filename, call. = FALSE)
  }

  data$include_primary_analysis <- figure_s1_primary_flag(
    data$include_primary_analysis,
    filename
  )
  data
}

prepare_figure_s1_alpha_data <- function(
    project_root = find_figure_s1_project_root()) {
  files <- figure_s1_source_files(project_root)
  all_files <- unlist(files, use.names = FALSE)
  if (!all(file.exists(all_files))) {
    stop(
      "Missing Figure S1 source file(s): ",
      paste(all_files[!file.exists(all_files)], collapse = ", "),
      call. = FALSE
    )
  }

  alpha_sources <- list(
    read_figure_s1_alpha_source(files$alpha[[1]], "SampleID"),
    read_figure_s1_alpha_source(files$alpha[[2]], "sampleID")
  )
  manifests <- Map(
    read_figure_s1_manifest,
    unname(files$metadata),
    unname(figure_s1_batch_ids)
  )
  names(alpha_sources) <- unname(figure_s1_batch_ids)
  names(manifests) <- unname(figure_s1_batch_ids)

  source_ids_without_manifest <- Map(function(alpha, manifest) {
    sort(setdiff(alpha$sample_id, manifest$sample_id))
  }, alpha_sources, manifests)
  excluded_nonprimary_ids <- lapply(manifests, function(manifest) {
    sort(manifest$sample_id[!manifest$include_primary_analysis])
  })

  prepared <- Map(function(alpha, manifest) {
    primary_manifest <- manifest[manifest$include_primary_analysis, , drop = FALSE]
    missing_alpha_ids <- setdiff(primary_manifest$sample_id, alpha$sample_id)
    if (length(missing_alpha_ids) > 0L) {
      stop(
        "Primary manifest samples are absent from an alpha-diversity source: ",
        paste(missing_alpha_ids, collapse = ", "),
        call. = FALSE
      )
    }

    alpha_rows <- alpha[
      match(primary_manifest$sample_id, alpha$sample_id),
      ,
      drop = FALSE
    ]
    metadata_rows <- primary_manifest[
      ,
      setdiff(names(primary_manifest), "sample_id"),
      drop = FALSE
    ]
    cbind(alpha_rows, metadata_rows)
  }, alpha_sources, manifests)

  result <- do.call(rbind, prepared)
  rownames(result) <- NULL
  result$treatment <- factor(
    result$treatment,
    levels = figure_s1_treatment_levels
  )
  result$sex <- factor(result$sex, levels = figure_s1_sex_levels)
  attr(result, "source_ids_without_manifest") <- source_ids_without_manifest
  attr(result, "excluded_nonprimary_ids") <- excluded_nonprimary_ids
  result
}

prepare_figure_s1_panel_data <- function(alpha_data) {
  cohort_data <- lapply(figure_s1_cohort_configurations, function(configuration) {
    selected <- alpha_data[
      alpha_data$batch_id == configuration$batch_id &
        alpha_data$generation == configuration$generation,
      ,
      drop = FALSE
    ]
    if (nrow(selected) != configuration$expected_sample_count) {
      stop(
        "Unexpected sample count for Figure S1 cohort ",
        configuration$label,
        ".",
        call. = FALSE
      )
    }
    selected
  })

  lapply(figure_s1_panel_configurations, function(configuration) {
    cohort_data[[configuration$cohort]]
  })
}

calculate_figure_s1_statistics <- function(data, metric, method) {
  if (!metric %in% c("Observed", "InvSimpson")) {
    stop("Unsupported Figure S1 metric: ", metric, call. = FALSE)
  }
  if (!method %in% c("parametric", "nonparametric")) {
    stop("Unsupported Figure S1 statistical method: ", method, call. = FALSE)
  }

  sex_formula <- stats::reformulate("sex", response = metric)
  treatment_formula <- stats::reformulate("treatment", response = metric)
  sex_test <- NULL
  if (dplyr::n_distinct(data$sex) > 1L) {
    sex_test <- if (method == "parametric") {
      rstatix::t_test(data, sex_formula)
    } else {
      rstatix::wilcox_test(data, sex_formula)
    }
  }

  grouped <- dplyr::group_by(data, sex, .drop = TRUE)
  treatment_test <- if (method == "parametric") {
    rstatix::tukey_hsd(grouped, treatment_formula)
  } else {
    rstatix::dunn_test(grouped, treatment_formula, p.adjust.method = "fdr")
  }

  significant_test <- dplyr::filter(treatment_test, p.adj.signif != "ns")
  if (nrow(significant_test) > 0L) {
    significant_test <- rstatix::add_y_position(significant_test)
  }

  list(
    sex = sex_test,
    treatment = treatment_test,
    significant = significant_test
  )
}

figure_s1_pwc_label <- function(treatment_test) {
  test_arguments <- attr(treatment_test, "args")
  if (identical(test_arguments$method, "tukey_hsd")) {
    return("pwc: Tukey HSD; p.adjust: Tukey")
  }
  if (identical(test_arguments$method, "dunn_test")) {
    return("pwc: Dunn test; p.adjust: Fdr")
  }
  stop("Unsupported Figure S1 pairwise-comparison result.", call. = FALSE)
}

make_figure_s1_panel <- function(data, metric, statistics) {
  subtitle <- NULL
  if (!is.null(statistics$sex)) {
    subtitle <- paste(
      "Sex comparison: ",
      rstatix::get_test_label(statistics$sex, detailed = TRUE, type = "text")
    )
  }

  plot <- ggpubr::ggboxplot(
    data,
    x = "treatment",
    y = metric,
    fill = "treatment",
    facet.by = "sex",
    palette = unname(figure_s1_colors),
    outlier.shape = NA
  ) +
    ggplot2::labs(
      title = metric,
      subtitle = subtitle,
      caption = figure_s1_pwc_label(statistics$treatment),
      x = NULL,
      y = metric,
      fill = "Group"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Group")) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.direction = "vertical",
      legend.position = "right",
      plot.subtitle = ggplot2::element_text(hjust = 1, size = 9),
      panel.background = ggplot2::element_rect(
        fill = "gray93",
        color = "grey",
        linetype = "solid"
      ),
      panel.grid.major = ggplot2::element_line(linewidth = 0.5, color = "white")
    )

  if (nrow(statistics$significant) > 0L) {
    plot <- plot + ggpubr::stat_pvalue_manual(
      statistics$significant,
      label = "p.adj.signif",
      tip.length = 0.01
    )
  }

  plot
}

generate_figure_s1_ah <- function(
    output_dir = file.path(
      find_figure_s1_project_root(),
      "supplementary_figures",
      "Figure_S1",
      "output"
    ),
    project_root = find_figure_s1_project_root()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  alpha_data <- prepare_figure_s1_alpha_data(project_root)
  panel_data <- prepare_figure_s1_panel_data(alpha_data)
  statistics <- Map(function(data, configuration) {
    calculate_figure_s1_statistics(
      data,
      configuration$metric,
      configuration$method
    )
  }, panel_data, figure_s1_panel_configurations)
  plots <- Map(function(data, configuration, panel_statistics) {
    make_figure_s1_panel(data, configuration$metric, panel_statistics)
  }, panel_data, figure_s1_panel_configurations, statistics)

  files <- setNames(
    file.path(
      output_dir,
      paste0("FigureS1", names(figure_s1_panel_configurations), "_generated.png")
    ),
    names(figure_s1_panel_configurations)
  )
  for (panel in names(files)) {
    cohort <- figure_s1_panel_configurations[[panel]]$cohort
    ggplot2::ggsave(
      filename = files[[panel]],
      plot = plots[[panel]],
      width = figure_s1_cohort_configurations[[cohort]]$width,
      height = 5,
      dpi = 150
    )
  }

  invisible(list(
    files = files,
    plots = plots,
    data = panel_data,
    statistics = statistics,
    configurations = figure_s1_panel_configurations
  ))
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s1_ah()
  message(
    "Generated Figure S1 panels A-H in: ",
    normalizePath(dirname(result$files[[1]]), winslash = "/")
  )
}
