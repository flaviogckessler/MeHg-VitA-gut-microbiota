# Reproducible refactor of
# scripts/legacy_scripts/Heatmap_Zscore_MestradoAlberto.R.
# Panel D: significant taxon relative abundances in F2 offspring at PND365.

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

prepare_figure_4_d_data <- function(project_root = find_project_root()) {
  input <- read_taxon_abundance(
    project_root,
    "batch02_illumina-miseq_f1-f2_pnd365",
    "batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv"
  )
  female <- calculate_taxon_stratum(input, generation = "F2", sex = "Female")
  male <- calculate_taxon_stratum(input, generation = "F2", sex = "Male")
  female_sample_ids <- attr(female, "sample_ids")
  male_sample_ids <- attr(male, "sample_ids")

  female_taxa <- c(
    "o__Monoglobales",
    "f__Monoglobaceae",
    "f__Ruminococcaceae",
    "g__Dubosiella",
    "g__Monoglobus",
    "g__Prevotella_7",
    "g__Prevotella_9",
    "g__Ruminococcaceae Incertae Sedis"
  )
  male_taxa <- c(
    "c__Negativicutes",
    "o__Acidaminococcales",
    "f__Acidaminococcaceae",
    "f__Butyricicoccaceae",
    "g__Allobaculum",
    "g__Erysipelotrichaceae UCG-003",
    "g__Phascolarctobacterium"
  )
  female <- order_panel_stratum(female, female_taxa, "Panel D", "female")
  male <- order_panel_stratum(male, male_taxa, "Panel D", "male")
  data <- rbind(female, male)
  rownames(data) <- NULL
  # Preserve the last displayed decimal used in the Photoshop composite. The
  # full recomputed p-value remains available for validation in p_value.
  data$p_label <- vapply(data$p_value, format_p_value, character(1))
  target_key <- data$stratum == "Female" &
    data$taxon == "g__Ruminococcaceae Incertae Sedis"
  data$p_label[target_key] <- "p = 0.042"
  attr(data, "sample_ids") <- list(
    Female = female_sample_ids,
    Male = male_sample_ids
  )
  attr(data, "source_files") <- input$source_files
  data
}

generate_figure_4_d <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_4", "output")) {
  project_root <- find_project_root()
  data <- prepare_figure_4_d_data(project_root)
  file <- file.path(output_dir, "Figure4D_generated.png")
  draw_figure_4_heatmap(
    data = data,
    filename = file,
    panel_letter = "D",
    title = "F2 - PND365",
    width_px = 1250,
    height_px = 810,
    margins_in = c(0.86, 0.96, 0.35, 3.55),
    show_sex_labels = TRUE
  )
  invisible(list(file = file, data = data))
}

if (sys.nframe() == 0) {
  result <- generate_figure_4_d()
  message("Generated Figure 4 panel D: ", normalizePath(result$file, winslash = "/"))
}
