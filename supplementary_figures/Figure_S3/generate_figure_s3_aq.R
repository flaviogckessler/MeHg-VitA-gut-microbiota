# Reproducible Figure S3 panels A-Q.
#
# Figure S3 contains the 10 female and seven male F1-PND21 taxon-abundance
# panels used in the supplied Photoshop composite. This script generates the
# individual panels only.

find_figure_s3_common_script <- function(start = getwd()) {
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

source(find_figure_s3_common_script())

figure_s3_labels <- c(
  "c__Alphaproteobacteria",
  "o__Peptococcales",
  "o__Rhodospirillales",
  "f__Peptococcaceae",
  "g__Dorea",
  "g__Flavonifractor",
  "g__Fusicatenibacter",
  "g__Lachnospiraceae.AC2044.group",
  "g__Lachnospiraceae.FCS020.group",
  "g__Prevotellaceae.NK3B31.group",
  "p__Cyanobacteria",
  "c__Vampirivibrionia",
  "o__Gastranaerophilales",
  "f__Aerococcaceae",
  "g__Aerococcus",
  "g__Lactococcus",
  "g__Prevotellaceae.NK3B31.group"
)

figure_s3_specification <- list(
  figure_id = "S3",
  batch_id = unname(supplementary_abundance_batch_ids[["batch1"]]),
  generation = "F1",
  panels = data.frame(
    panel = supplementary_abundance_panel_codes(length(figure_s3_labels)),
    sex = c(rep("Female", 10L), rep("Male", 7L)),
    label = figure_s3_labels,
    stringsAsFactors = FALSE
  ),
  output_subdirectories = c(Female = "females", Male = "male"),
  # The historical PND21 source panels plot the character-valued legacy Grupo
  # field in this order. The order and positional colors are preserved here so
  # the individual panels match the supplied reference exports exactly.
  treatment_levels = c("CN", "MeHg", "VitA", "VitA+MeHg"),
  treatment_colors = c(
    CN = "#98F5FF",
    MeHg = "#B3EE3A",
    VitA = "#CD2626",
    `VitA+MeHg` = "#4D4D4D"
  ),
  expected_sample_counts = c(Female = 18L, Male = 16L),
  expected_treatment_counts = list(
    Female = c(CN = 5L, MeHg = 5L, VitA = 4L, `VitA+MeHg` = 4L),
    Male = c(CN = 5L, MeHg = 5L, VitA = 3L, `VitA+MeHg` = 3L)
  ),
  composite_filename = "Figure_S3_FinalFigure.png",
  composite_dimensions = c(2790, 3169),
  expected_abundance_md5 = "92beb6c9b419284b9dc8cd3574cd7a86",
  expected_manifest_md5 = "17734fb98113e31be5d2b8834e832b24",
  expected_omnibus_statistics = c(
    7.98477970627504, 8.34912280701755, 7.98477970627504,
    8.34912280701755, 11.2459322033898, 9.29893992932862,
    10.7473684210526, 8.09658031088082, 10.0101159114858,
    7.85087719298245, 9.8060606060606, 9.8060606060606,
    9.8060606060606, 9.24444444444445, 9.24444444444445,
    9.24444444444445, 9.41764705882353
  ),
  expected_omnibus_p = c(
    0.0463273128101511, 0.0393214117666461, 0.0463273128101511,
    0.0393214117666461, 0.0104677136083723, 0.0255693625558909,
    0.0131734092849804, 0.044057287322268, 0.0184803419353083,
    0.0491962864921898, 0.0202887118600479, 0.0202887118600479,
    0.0202887118600479, 0.0262113670860976, 0.0262113670860976,
    0.0262113670860976, 0.0242237904287946
  ),
  expected_significant_counts = c(
    0L, 1L, 0L, 1L, 2L, 2L, 2L, 0L, 1L, 0L, 1L, 1L, 1L, 3L, 3L,
    3L, 1L
  ),
  expected_significant_signatures = c(
    "", "CN|VitA", "", "CN|VitA", "CN|MeHg,MeHg|VitA",
    "CN|MeHg,CN|VitA+MeHg", "CN|VitA,VitA|VitA+MeHg", "",
    "CN|VitA+MeHg", "", "MeHg|VitA", "MeHg|VitA", "MeHg|VitA",
    "CN|VitA,MeHg|VitA,VitA|VitA+MeHg",
    "CN|VitA,MeHg|VitA,VitA|VitA+MeHg",
    "CN|VitA,MeHg|VitA,VitA|VitA+MeHg", "MeHg|VitA"
  )
)

generate_figure_s3_aq <- function(
    output_dir = file.path(
      find_supplementary_abundance_project_root(),
      "supplementary_figures",
      "Figure_S3",
      "output"
    )) {
  generate_supplementary_abundance_figure(
    figure_s3_specification,
    output_dir
  )
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s3_aq()
  message(
    "Generated Figure S3 panels A-Q in: ",
    normalizePath(dirname(result$files[[1L]]), winslash = "/")
  )
}
