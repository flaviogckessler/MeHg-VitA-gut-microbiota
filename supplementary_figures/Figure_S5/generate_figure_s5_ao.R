# Reproducible Figure S5 panels A-O.
#
# Figure S5 contains the eight female and seven male F2-PND365 taxon-abundance
# panels used in the supplied Photoshop composite. This script generates the
# individual panels only.

# Load Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
})

find_figure_s5_common_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "figure_s2_s5_abundance_common.R"),
      file.path(
        path,
        "supplementary_figures",
        "figure_s2_s5_abundance_common.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate figure_s2_s5_abundance_common.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s5_common_script())

figure_s5_labels <- c(
  "o__Monoglobales",
  "f__Monoglobaceae",
  "f__Ruminococcaceae",
  "g__Dubosiella",
  "g__Monoglobus",
  "g__Prevotella_7",
  "g__Prevotella_9",
  "g__Ruminococcaceae.Incertae.Sedis",
  "c__Negativicutes",
  "o__Acidaminococcales",
  "f__Acidaminococcaceae",
  "f__Butyricicoccaceae",
  "g__Allobaculum",
  "g__Erysipelotrichaceae.UCG.003",
  "g__Phascolarctobacterium"
)

figure_s5_specification <- list(
  figure_id = "S5",
  batch_id = unname(supplementary_abundance_batch_ids[["batch2"]]),
  generation = "F2",
  panels = data.frame(
    panel = supplementary_abundance_panel_codes(length(figure_s5_labels)),
    sex = c(rep("Female", 8L), rep("Male", 7L)),
    label = figure_s5_labels,
    stringsAsFactors = FALSE
  ),
  output_subdirectories = c(Female = "females", Male = "male"),
  treatment_levels = c("CN", "VitA", "MeHg", "VitA+MeHg"),
  treatment_colors = c(
    CN = "#98F5FF",
    VitA = "#B3EE3A",
    MeHg = "#CD2626",
    `VitA+MeHg` = "#4D4D4D"
  ),
  expected_sample_counts = c(Female = 12L, Male = 12L),
  expected_treatment_counts = list(
    Female = c(CN = 3L, VitA = 3L, MeHg = 3L, `VitA+MeHg` = 3L),
    Male = c(CN = 3L, VitA = 3L, MeHg = 3L, `VitA+MeHg` = 3L)
  ),
  composite_filename = "Figure_S5_FinalFigure.png",
  composite_dimensions = c(4668, 4044),
  expected_abundance_md5 = "f3f65293ad4e35e647e6d93d393c5afb",
  expected_manifest_md5 = "a703799a2d60b1c3d935901ad76d8d6c",
  expected_omnibus_statistics = c(
    8.94871794871795, 8.94871794871795, 8.07692307692308,
    8.18796992481203, 8.94871794871795, 8.29537953795379,
    8.29537953795379, 8.23076923076923, 8.19420289855073,
    8.19420289855073, 8.19420289855073, 9.46153846153846,
    8.12820512820513, 8.60424966799469, 8.19420289855073
  ),
  expected_omnibus_p = c(
    0.0299805369160345, 0.0299805369160345, 0.044448402794751,
    0.042282545635888, 0.0299805369160345, 0.0402857078711178,
    0.0402857078711178, 0.0414755598461637, 0.0421640767150978,
    0.0421640767150978, 0.0421640767150978, 0.0237440664767199,
    0.0434350841537281, 0.0350427187131641, 0.0421640767150978
  ),
  expected_significant_counts = c(
    0L, 0L, 1L, 1L, 0L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L
  ),
  expected_significant_signatures = c(
    "", "", "MeHg|VitA+MeHg", "VitA|VitA+MeHg", "",
    "CN|VitA+MeHg,VitA|VitA+MeHg", "CN|VitA+MeHg,VitA|VitA+MeHg",
    "MeHg|VitA+MeHg", "VitA|MeHg,MeHg|VitA+MeHg",
    "VitA|MeHg,MeHg|VitA+MeHg", "VitA|MeHg,MeHg|VitA+MeHg",
    "CN|VitA", "CN|VitA", "VitA|VitA+MeHg",
    "VitA|MeHg,MeHg|VitA+MeHg"
  )
)

generate_figure_s5_ao <- function(
    output_dir = file.path(
      find_supplementary_abundance_project_root(),
      "supplementary_figures",
      "Figure_S5",
      "output"
    )) {
  generate_supplementary_abundance_figure(
    figure_s5_specification,
    output_dir
  )
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s5_ao()
  message(
    "Generated Figure S5 panels A-O in: ",
    normalizePath(dirname(result$files[[1L]]), winslash = "/")
  )
}
