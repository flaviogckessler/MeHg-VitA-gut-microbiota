# Reproducible refactor of
# scripts/legacy_scripts/Abundancias_Barras_Mestrado.R.
# The legacy file is retained for provenance and is not a runtime dependency.

# - Packages
library(tidyverse)

# - Repository paths
find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    markers <- c(
      file.path(path, "data", "processed"),
      file.path(path, "data", "metadata"),
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

# - Panel A data
prepare_panel_a_data <- function(project_root) {
  count_file <- file.path(
    project_root,
    "data", "processed",
    "batch01_iontorrent_dams-f0_f1-pnd21",
    "asv_counts_with_taxonomy.tsv"
  )
  metadata_file <- file.path(
    project_root,
    "data", "metadata",
    "batch01_iontorrent_dams-f0_f1-pnd21",
    "sample_manifest.csv"
  )

  counts <- read_tsv(count_file, show_col_types = FALSE)
  metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
    filter(include_primary_analysis)
  sample_ids <- intersect(names(counts), metadata$sample_id)

  data <- counts %>%
    select(OTU, Phylum, all_of(sample_ids)) %>%
    pivot_longer(
      cols = all_of(sample_ids),
      names_to = "sample_id",
      values_to = "count"
    ) %>%
    left_join(
      metadata %>% select(sample_id, generation, treatment, sex),
      by = "sample_id"
    ) %>%
    mutate(
      generation_sex = case_when(
        generation == "F0" ~ "F0_Dams",
        TRUE ~ paste(generation, sex, sep = "_")
      )
    )

  phylum_levels <- data %>%
    group_by(sample_id, Phylum) %>%
    summarise(phylum_count = sum(count), .groups = "drop") %>%
    group_by(sample_id) %>%
    mutate(phylum_frequency = phylum_count / sum(phylum_count)) %>%
    ungroup() %>%
    group_by(Phylum) %>%
    summarise(mean_frequency = mean(phylum_frequency), .groups = "drop") %>%
    arrange(desc(mean_frequency)) %>%
    pull(Phylum)

  data %>%
    mutate(
      Phylum = factor(Phylum, levels = phylum_levels[!is.na(phylum_levels)]),
      treatment = factor(treatment, levels = c("CN", "VitA", "MeHg", "VitA+MeHg")),
      generation_sex = factor(
        generation_sex,
        levels = c("F0_Dams", "F1_Female", "F1_Male")
      )
    )
}

# - Panels B and C data
prepare_panels_bc_data <- function(project_root) {
  abundance_file <- file.path(
    project_root,
    "data", "processed",
    "batch02_illumina-miseq_f1-f2_pnd365",
    "asv_relative_abundance_with_taxonomy.tsv"
  )
  metadata_file <- file.path(
    project_root,
    "data", "metadata",
    "batch02_illumina-miseq_f1-f2_pnd365",
    "sample_manifest.csv"
  )

  abundance <- read_tsv(abundance_file, show_col_types = FALSE)
  metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
    filter(include_primary_analysis)
  sample_ids <- intersect(names(abundance), metadata$sample_id)
  phylum_levels <- c(
    "Firmicutes", "Bacteroidota", "Proteobacteria", "Actinobacteriota",
    "Campylobacterota", "Desulfobacterota", "Verrucomicrobiota",
    "Patescibacteria"
  )

  abundance %>%
    select(ASV, Phylum, all_of(sample_ids)) %>%
    pivot_longer(
      cols = all_of(sample_ids),
      names_to = "sample_id",
      values_to = "abundance"
    ) %>%
    group_by(sample_id, Phylum) %>%
    summarise(abundance = sum(abundance), .groups = "drop") %>%
    left_join(
      metadata %>% select(sample_id, generation, treatment, sex),
      by = "sample_id"
    ) %>%
    mutate(
      Phylum = factor(Phylum, levels = phylum_levels),
      treatment = factor(treatment, levels = c("CN", "VitA", "MeHg", "VitA+MeHg")),
      sex = factor(sex, levels = c("Female", "Male"))
    )
}

# - Panel plots
make_panel_a <- function(data) {
  ggplot(data, aes(x = sample_id, y = count)) +
    facet_grid(
      cols = vars(generation_sex, treatment),
      scales = "free_x",
      space = "free_x"
    ) +
    geom_col(aes(fill = Phylum), position = "fill", width = 1) +
    scale_fill_brewer(palette = "Paired") +
    scale_y_continuous(
      name = "Relative abundance",
      labels = scales::label_percent()
    ) +
    labs(x = "SampleID", fill = "Phylum") +
    theme(
      axis.text.x = element_text(angle = 90, size = 6),
      axis.text.y = element_text(color = "black"),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.clip = "off"
    )
}

make_panel_bc <- function(data, selected_generation) {
  data %>%
    filter(generation == selected_generation) %>%
    ggplot(aes(x = sample_id, y = abundance, fill = Phylum)) +
    geom_col(width = 1) +
    facet_wrap(
      vars(sex, treatment),
      scales = "free_x",
      nrow = 1
    ) +
    scale_fill_brewer(palette = "Paired", drop = FALSE) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(
      title = "2739 Taxa (Phylum)",
      x = "SampleID",
      y = "Relative abundance",
      fill = "Phylum"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, size = 6),
      axis.text.y = element_text(color = "black"),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank()
    )
}

# - Figure generation
generate_figure_3_abc <- function(
    output_dir = file.path(
      find_project_root(), "figures", "Figure_3", "output"
    )) {
  project_root <- find_project_root()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  panel_a_data <- prepare_panel_a_data(project_root)
  panels_bc_data <- prepare_panels_bc_data(project_root)
  plots <- list(
    A = make_panel_a(panel_a_data),
    B = make_panel_bc(panels_bc_data, "F1"),
    C = make_panel_bc(panels_bc_data, "F2")
  )
  files <- c(
    A = file.path(output_dir, "Figure3A_generated.png"),
    B = file.path(output_dir, "Figure3B_generated.png"),
    C = file.path(output_dir, "Figure3C_generated.png")
  )

  ggsave(files[["A"]], plots[["A"]], width = 12, height = 7, dpi = 150)
  ggsave(files[["B"]], plots[["B"]], width = 10, height = 8, dpi = 150)
  ggsave(files[["C"]], plots[["C"]], width = 10, height = 8, dpi = 150)

  invisible(list(files = files, plots = plots))
}

if (sys.nframe() == 0) {
  result <- generate_figure_3_abc()
  message("Generated Figure 3 panels A-C in: ", dirname(result$files[[1]]))
}
