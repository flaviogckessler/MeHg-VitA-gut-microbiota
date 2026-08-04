# - Packages
library(tidyverse)
library(rstatix)

# - Repository paths
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "metadata"),
      file.path(path, "results", "taxon_relative_abundance"),
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

# - Canonical abundance data
phylum_levels <- c(
  "Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota"
)
treatment_levels <- c("CN", "VitA", "MeHg", "VitA+MeHg")
treatment_colors <- c(
  CN = "#98F5FF",
  VitA = "#B3EE3A",
  MeHg = "#CD2626",
  `VitA+MeHg` = "#4D4D4D"
)

read_batch_abundance <- function(project_root, batch_id, filename) {
  abundance <- read_csv(
    file.path(project_root, "results", "taxon_relative_abundance", filename),
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  expected_taxa <- paste0("p__", phylum_levels)
  phylum_rows <- abundance %>%
    filter(regex_match %in% expected_taxa)

  if (nrow(phylum_rows) != length(expected_taxa) ||
      !setequal(phylum_rows$regex_match, expected_taxa)) {
    stop("The abundance table must contain one row for each required phylum.",
         call. = FALSE)
  }

  descriptor_columns <- c(
    "Num", "taxon_id", "tax_rank", "tax_name", "regex_match"
  )
  sample_columns <- setdiff(names(abundance), descriptor_columns)
  abundance_long <- phylum_rows %>%
    transmute(
      phylum = str_remove(regex_match, "^p__"),
      across(all_of(sample_columns))
    ) %>%
    pivot_longer(
      cols = all_of(sample_columns),
      names_to = "sample_id",
      values_to = "relative_abundance"
    )

  metadata <- read_csv(
    file.path(project_root, "data", "metadata", batch_id, "sample_manifest.csv"),
    show_col_types = FALSE
  ) %>%
    filter(include_primary_analysis) %>%
    select(
      sample_id, batch_id, generation, life_stage, collection_day_pnd,
      treatment, sex
    )

  missing_samples <- setdiff(metadata$sample_id, abundance_long$sample_id)
  if (length(missing_samples) > 0) {
    stop(
      paste("Included samples missing from the abundance table:",
            paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  inner_join(abundance_long, metadata, by = "sample_id") %>%
    mutate(
      phylum = factor(phylum, levels = phylum_levels),
      treatment = factor(treatment, levels = treatment_levels),
      sex = factor(sex, levels = c("Female", "Male"))
    )
}

prepare_abundance_data <- function(project_root) {
  list(
    batch1 = read_batch_abundance(
      project_root,
      "batch01_iontorrent_dams-f0_f1-pnd21",
      "batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv"
    ),
    batch2 = read_batch_abundance(
      project_root,
      "batch02_illumina-miseq_f1-f2_pnd365",
      "batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv"
    )
  )
}

# - Panel subsets
prepare_panel_hk_data <- function(abundance_data) {
  list(
    H = abundance_data$batch1 %>% filter(generation == "F0"),
    I = abundance_data$batch1 %>% filter(generation == "F1"),
    J = abundance_data$batch2 %>% filter(generation == "F1"),
    K = abundance_data$batch2 %>% filter(generation == "F2")
  )
}

# - Statistical comparisons
calculate_panel_hk_statistics <- function(data) {
  shapiro <- data %>%
    group_by(phylum) %>%
    shapiro_test(relative_abundance) %>%
    ungroup() %>%
    mutate(
      method = if_else(p > 0.05, "parametric", "nonparametric")
    )

  treatment <- map_dfr(phylum_levels, function(current_phylum) {
    phylum_data <- data %>% filter(phylum == current_phylum)
    method <- shapiro %>%
      filter(phylum == current_phylum) %>%
      pull(method)

    result <- if (method == "parametric") {
      tukey_hsd(phylum_data, relative_abundance ~ treatment)
    } else {
      dunn_test(
        phylum_data,
        relative_abundance ~ treatment,
        p.adjust.method = "fdr"
      )
    }

    result %>%
      mutate(
        phylum = factor(current_phylum, levels = phylum_levels),
        method = method,
        .before = 1
      )
  })

  sex <- NULL
  if (n_distinct(data$sex) > 1) {
    sex <- data %>%
      group_by(phylum) %>%
      wilcox_test(relative_abundance ~ sex) %>%
      ungroup()
  }

  list(shapiro = shapiro, treatment = treatment, sex = sex)
}

# - Plot labels
join_names <- function(values) {
  if (length(values) == 1) return(values)
  paste0(paste(values[-length(values)], collapse = ", "),
         " and ", values[[length(values)]])
}

make_shapiro_caption <- function(shapiro) {
  normal <- as.character(shapiro$phylum[shapiro$p > 0.05])
  nonnormal <- as.character(shapiro$phylum[shapiro$p <= 0.05])

  if (length(normal) == 0) {
    return("Shapiro-Wilk tests: all phyla p < 0.05")
  }
  if (length(nonnormal) == 0) {
    return("Shapiro-Wilk tests: all phyla p > 0.05")
  }

  paste0(
    "Shapiro-Wilk tests: ", join_names(normal), " p > 0.05; ",
    join_names(nonnormal), " p < 0.05"
  )
}

make_sex_subtitle <- function(sex_tests) {
  if (is.null(sex_tests)) return(NULL)
  if (all(sex_tests$p > 0.05)) {
    return("Sex comparisons for all phyla: Wilcoxon tests p > 0.05")
  }

  significant <- as.character(sex_tests$phylum[sex_tests$p <= 0.05])
  paste0("Sex comparisons: Wilcoxon tests p <= 0.05 for ",
         join_names(significant))
}

# - Panel plots
make_phylum_panel <- function(data, title) {
  statistics <- calculate_panel_hk_statistics(data)

  plot <- ggplot(
    data,
    aes(x = treatment, y = relative_abundance, fill = treatment)
  ) +
    geom_boxplot(width = 0.72) +
    scale_fill_manual(values = treatment_colors, drop = FALSE) +
    facet_wrap(~ phylum, scales = "free_y", nrow = 2) +
    labs(
      title = title,
      subtitle = make_sex_subtitle(statistics$sex),
      caption = str_wrap(
        make_shapiro_caption(statistics$shapiro),
        width = 78
      ),
      x = "Phylum",
      y = "Relative Abundance",
      fill = "Group"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.direction = "vertical",
      legend.position = "right",
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(
        size = 8, hjust = 0.5, lineheight = 1, margin = margin(t = 6)
      ),
      plot.margin = margin(6, 8, 10, 8),
      panel.background = element_rect(
        fill = "gray93", color = "grey", linetype = "solid"
      ),
      panel.grid.major = element_line(linewidth = 0.5, color = "white")
    )

  list(plot = plot, statistics = statistics)
}

# - Figure generation
generate_figure_3_hk <- function(
    output_dir = file.path(
      find_project_root(), "figures", "Figure_3", "output"
    )) {
  project_root <- find_project_root()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  panel_data <- prepare_panel_hk_data(prepare_abundance_data(project_root))
  titles <- c(
    H = "Dams (F0)",
    I = "F1 - PND21",
    J = "F1 - PND365",
    K = "F2 - PND365"
  )
  panels <- map2(panel_data, titles, make_phylum_panel)
  plots <- map(panels, "plot")
  statistics <- map(panels, "statistics")
  files <- set_names(
    file.path(output_dir, paste0("Figure3", names(panel_data), "_generated.png")),
    names(panel_data)
  )

  walk2(files, plots, function(filename, plot) {
    ggsave(filename, plot, width = 6, height = 5, dpi = 150)
  })

  invisible(list(
    files = files,
    plots = plots,
    data = panel_data,
    statistics = statistics
  ))
}

if (sys.nframe() == 0) {
  result <- generate_figure_3_hk()
  message(
    "Generated Figure 3 panels H-K in: ",
    normalizePath(dirname(result$files[[1]]), winslash = "/")
  )
}
