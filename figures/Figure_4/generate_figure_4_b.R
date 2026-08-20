# Reproducible refactor of
# scripts/legacy_scripts/Heatmap_Zscore_MestradoAlberto.R.
# Panel B: significant taxon relative abundances in F1 offspring at PND21.
# Shared base-R input, statistical, and heatmap helpers are kept in the copied
# panel-A workflow so that Figure 4 has exactly one generation script per panel.

# Load Packages ----
suppressPackageStartupMessages({
  library(graphics)
  library(grDevices)
  library(stats)
})

find_panel_a_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "generate_figure_4_a.R"),
      file.path(path, "figures", "Figure_4", "generate_figure_4_a.R")
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_4_a.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_panel_a_script(), local = TRUE)

order_panel_stratum <- function(data, expected_taxa, panel, stratum) {
  if (!setequal(data$taxon, expected_taxa)) {
    stop(
      paste(panel, stratum, "significant taxa differ from the final-figure specification."),
      call. = FALSE
    )
  }
  data[match(expected_taxa, data$taxon), , drop = FALSE]
}

prepare_figure_4_b_data <- function(project_root = find_project_root()) {
  input <- read_taxon_abundance(
    project_root,
    "batch01_iontorrent_dams-f0_f1-pnd21",
    "batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv"
  )
  female <- calculate_taxon_stratum(input, generation = "F1", sex = "Female")
  male <- calculate_taxon_stratum(input, generation = "F1", sex = "Male")
  female_sample_ids <- attr(female, "sample_ids")
  male_sample_ids <- attr(male, "sample_ids")

  female_taxa <- c(
    "c__Alphaproteobacteria",
    "o__Peptococcales",
    "o__Rhodospirillales",
    "f__Peptococcaceae",
    "g__Dorea",
    "g__Flavonifractor",
    "g__Fusicatenibacter",
    "g__Lachnospiraceae AC2044 group",
    "g__Lachnospiraceae FCS020 group",
    "g__Prevotellaceae NK3B31 group"
  )
  male_taxa <- c(
    "p__Cyanobacteria",
    "c__Vampirivibrionia",
    "o__Gastranaerophilales",
    "f__Aerococcaceae",
    "g__Aerococcus",
    "g__Lactococcus",
    "g__Prevotellaceae NK3B31 group"
  )
  female <- order_panel_stratum(female, female_taxa, "Panel B", "female")
  male <- order_panel_stratum(male, male_taxa, "Panel B", "male")
  data <- rbind(female, male)
  rownames(data) <- NULL
  attr(data, "sample_ids") <- list(
    Female = female_sample_ids,
    Male = male_sample_ids
  )
  attr(data, "source_files") <- input$source_files
  data
}

generate_figure_4_b <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_4", "output")) {
  project_root <- find_project_root()
  data <- prepare_figure_4_b_data(project_root)
  file <- file.path(output_dir, "Figure4B_generated.png")
  draw_figure_4_heatmap(
    data = data,
    filename = file,
    panel_letter = "B",
    title = "F1 - PND21",
    width_px = 1150,
    height_px = 825,
    margins_in = c(0.86, 0.96, 0.35, 3.005),
    show_sex_labels = TRUE
  )
  invisible(list(file = file, data = data))
}

if (sys.nframe() == 0) {
  result <- generate_figure_4_b()
  message("Generated Figure 4 panel B: ", normalizePath(result$file, winslash = "/"))
}
