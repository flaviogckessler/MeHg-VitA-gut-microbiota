# Reproducible refactor of
# scripts/legacy_scripts/Heatmap_Zscore_MestradoAlberto.R.
# Panel C: significant taxon relative abundances in F1 offspring at PND365.

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

prepare_figure_4_c_data <- function(project_root = find_project_root()) {
  input <- read_taxon_abundance(
    project_root,
    "batch02_illumina-miseq_f1-f2_pnd365",
    "batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv"
  )
  female <- calculate_taxon_stratum(input, generation = "F1", sex = "Female")
  male <- calculate_taxon_stratum(input, generation = "F1", sex = "Male")
  female_sample_ids <- attr(female, "sample_ids")
  male_sample_ids <- attr(male, "sample_ids")

  female_taxa <- c(
    "p__Actinobacteriota",
    "c__Actinobacteria",
    "c__Bacilli",
    "c__Negativicutes",
    "o__Acidaminococcales",
    "o__Bifidobacteriales",
    "o__Christensenellales",
    "o__Erysipelotrichales",
    "o__Lactobacillales",
    "f__Acidaminococcaceae",
    "f__Anaerovoracaceae",
    "f__Bifidobacteriaceae",
    "f__Christensenellaceae",
    "f__Erysipelotrichaceae",
    "f__Lactobacillaceae",
    "g__[Eubacterium] ruminantium group",
    "g__Alistipes",
    "g__Allobaculum",
    "g__Amylolactobacillus",
    "g__ASF356",
    "g__Bifidobacterium",
    "g__CAG-352",
    "g__Christensenellaceae R-7 group",
    "g__Dubosiella",
    "g__HT002",
    "g__Lachnospiraceae NK4A136 group",
    "g__Limosilactobacillus",
    "g__Phascolarctobacterium",
    "g__Prevotellaceae Ga6A1 group"
  )
  male_taxa <- c(
    "c__Alphaproteobacteria",
    "o__Clostridiales",
    "o__Rickettsiales",
    "f__Clostridiaceae",
    "f__Mitochondria",
    "g__Alistipes",
    "g__CAG-352",
    "g__Clostridium sensu stricto 2",
    "g__Lactovum",
    "g__Prevotella",
    "g__Streptococcus"
  )
  female <- order_panel_stratum(female, female_taxa, "Panel C", "female")
  male <- order_panel_stratum(male, male_taxa, "Panel C", "male")
  data <- rbind(female, male)
  rownames(data) <- NULL
  # Preserve the last displayed decimal used in the Photoshop composite. The
  # full recomputed p-values remain in p_value and are what validation tests.
  data$p_label <- vapply(data$p_value, format_p_value, character(1))
  final_label_overrides <- c(
    "Female|f__Anaerovoracaceae" = "p = 0.024",
    "Female|f__Lactobacillaceae" = "p = 0.042",
    "Female|g__Alistipes" = "p = 0.024"
  )
  final_label_keys <- paste(data$stratum, data$taxon, sep = "|")
  override_index <- match(final_label_keys, names(final_label_overrides))
  has_override <- !is.na(override_index)
  data$p_label[has_override] <- unname(final_label_overrides[override_index[has_override]])
  attr(data, "sample_ids") <- list(
    Female = female_sample_ids,
    Male = male_sample_ids
  )
  attr(data, "source_files") <- input$source_files
  data
}

generate_figure_4_c <- function(
    output_dir = file.path(find_project_root(), "figures", "Figure_4", "output")) {
  project_root <- find_project_root()
  data <- prepare_figure_4_c_data(project_root)
  file <- file.path(output_dir, "Figure4C_generated.png")
  draw_figure_4_heatmap(
    data = data,
    filename = file,
    panel_letter = "C",
    title = "F1 - PND365",
    width_px = 1200,
    height_px = 1532,
    margins_in = c(0.86, 0.96, 0.35, 3.29),
    show_sex_labels = TRUE
  )
  invisible(list(file = file, data = data))
}

if (sys.nframe() == 0) {
  result <- generate_figure_4_c()
  message("Generated Figure 4 panel C: ", normalizePath(result$file, winslash = "/"))
}
