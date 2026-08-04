# - Packages
library(tidyverse)
library(rstatix)
library(ggpubr)

# - Repository paths
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata"),
      file.path(path, "results", "alpha_diversity"),
      file.path(path, "figures", "Figure_3")
    )
    if (all(dir.exists(markers))) return(path)

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the repository root.", call. = FALSE)
    }
    path <- parent
  }
}

# - Canonical alpha-diversity data
prepare_alpha_data <- function(project_root) {
  batch1_alpha <- read_csv(
    file.path(
      project_root, "results", "alpha_diversity",
      "batch01_iontorrent_dams-f0_f1-pnd21_alpha_diversity_with_metadata.csv"
    ),
    show_col_types = FALSE,
    name_repair = "unique_quiet"
  ) %>%
    transmute(sample_id = SampleID, Shannon)

  batch2_alpha <- read_csv(
    file.path(
      project_root, "results", "alpha_diversity",
      "batch02_illumina-miseq_f1-f2_pnd365_alpha_diversity_with_metadata.csv"
    ),
    show_col_types = FALSE
  ) %>%
    transmute(sample_id = sampleID, Shannon)

  batch1_metadata <- read_csv(
    file.path(
      project_root, "data", "metadata",
      "batch01_iontorrent_dams-f0_f1-pnd21", "sample_manifest.csv"
    ),
    show_col_types = FALSE
  ) %>%
    filter(include_primary_analysis)

  batch2_metadata <- read_csv(
    file.path(
      project_root, "data", "metadata",
      "batch02_illumina-miseq_f1-f2_pnd365", "sample_manifest.csv"
    ),
    show_col_types = FALSE
  ) %>%
    filter(include_primary_analysis)

  bind_rows(
    inner_join(batch1_alpha, batch1_metadata, by = "sample_id"),
    inner_join(batch2_alpha, batch2_metadata, by = "sample_id")
  ) %>%
    mutate(
      treatment = factor(
        treatment,
        levels = c("CN", "VitA", "MeHg", "VitA+MeHg")
      ),
      sex = factor(sex, levels = c("Female", "Male"))
    )
}

# - Panel subsets
prepare_panel_dg_data <- function(alpha_data) {
  list(
    D = alpha_data %>%
      filter(
        batch_id == "batch01_iontorrent_dams-f0_f1-pnd21",
        generation == "F0"
      ),
    E = alpha_data %>%
      filter(
        batch_id == "batch01_iontorrent_dams-f0_f1-pnd21",
        generation == "F1"
      ),
    F = alpha_data %>%
      filter(
        batch_id == "batch02_illumina-miseq_f1-f2_pnd365",
        generation == "F1"
      ),
    G = alpha_data %>%
      filter(
        batch_id == "batch02_illumina-miseq_f1-f2_pnd365",
        generation == "F2"
      )
  )
}

# - Statistical comparisons
calculate_panel_statistics <- function(data, method) {
  sex_test <- NULL
  if (n_distinct(data$sex) > 1) {
    sex_test <- if (method == "parametric") {
      t_test(data, Shannon ~ sex)
    } else {
      wilcox_test(data, Shannon ~ sex)
    }
  }

  treatment_test <- if (method == "parametric") {
    data %>%
      group_by(sex) %>%
      tukey_hsd(Shannon ~ treatment)
  } else {
    data %>%
      group_by(sex) %>%
      dunn_test(Shannon ~ treatment, p.adjust.method = "fdr")
  }

  significant_test <- treatment_test %>%
    filter(p.adj.signif != "ns")
  if (nrow(significant_test) > 0) {
    significant_test <- add_y_position(significant_test)
  }

  list(
    sex = sex_test,
    treatment = treatment_test,
    significant = significant_test
  )
}

# - Panel plots
make_shannon_panel <- function(data, method) {
  colors <- c(
    CN = "#98F5FF",
    VitA = "#B3EE3A",
    MeHg = "#CD2626",
    `VitA+MeHg` = "#4D4D4D"
  )
  statistics <- calculate_panel_statistics(data, method)
  subtitle <- NULL
  if (!is.null(statistics$sex)) {
    subtitle <- paste(
      "Sex comparison: ",
      get_test_label(statistics$sex, detailed = TRUE, type = "text")
    )
  }

  plot <- ggboxplot(
    data,
    x = "treatment",
    y = "Shannon",
    fill = "treatment",
    facet.by = "sex",
    palette = unname(colors),
    outlier.shape = NA
  ) +
    labs(
      title = "Shannon",
      subtitle = subtitle,
      caption = get_pwc_label(statistics$treatment),
      x = NULL,
      y = "Shannon",
      fill = "Group"
    ) +
    guides(fill = guide_legend(title = "Group")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.direction = "vertical",
      legend.position = "right",
      plot.subtitle = element_text(hjust = 1, size = 9),
      panel.background = element_rect(
        fill = "gray93",
        color = "grey",
        linetype = "solid"
      ),
      panel.grid.major = element_line(linewidth = 0.5, color = "white")
    )

  if (nrow(statistics$significant) > 0) {
    plot <- plot + stat_pvalue_manual(
      statistics$significant,
      label = "p.adj.signif",
      tip.length = 0.01
    )
  }

  plot
}

# - Figure generation
generate_figure_3_dg <- function(
    output_dir = file.path(
      find_project_root(), "figures", "Figure_3", "output"
    )) {
  project_root <- find_project_root()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  panel_data <- prepare_panel_dg_data(prepare_alpha_data(project_root))
  methods <- c(
    D = "parametric",
    E = "nonparametric",
    F = "nonparametric",
    G = "parametric"
  )
  plots <- map2(panel_data, methods, make_shannon_panel)
  files <- c(
    D = file.path(output_dir, "Figure3D_generated.png"),
    E = file.path(output_dir, "Figure3E_generated.png"),
    F = file.path(output_dir, "Figure3F_generated.png"),
    G = file.path(output_dir, "Figure3G_generated.png")
  )

  ggsave(files[["D"]], plots[["D"]], width = 4, height = 5, dpi = 150)
  walk(c("E", "F", "G"), function(panel) {
    ggsave(files[[panel]], plots[[panel]], width = 8, height = 5, dpi = 150)
  })

  invisible(list(files = files, plots = plots, data = panel_data))
}

if (sys.nframe() == 0) {
  result <- generate_figure_3_dg()
  message("Generated Figure 3 panels D-G in: ", dirname(result$files[[1]]))
}
