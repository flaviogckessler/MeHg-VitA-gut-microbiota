# Reproducible Figure S2 panels A-C.
#
# The supplied Figure_S2_FinalFigure.png was assembled manually in Adobe
# Photoshop. This script uses the shared Figures S2-S5 abundance workflow to
# generate only the three individual Dams (F0) relative-abundance boxplots.

find_figure_s2_common_script <- function(start = getwd()) {
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

source(find_figure_s2_common_script())

figure_s2_labels <- c(
  "g__.Eubacterium..ventriosum.group",
  "g__.Ruminococcus..gauvreauii.group",
  "g__Oribacterium"
)

figure_s2_specification <- list(
  figure_id = "S2",
  batch_id = unname(supplementary_abundance_batch_ids[["batch1"]]),
  generation = "F0",
  panels = data.frame(
    panel = supplementary_abundance_panel_codes(length(figure_s2_labels)),
    sex = rep("Female", length(figure_s2_labels)),
    label = figure_s2_labels,
    stringsAsFactors = FALSE
  ),
  treatment_levels = c("CN", "VitA", "MeHg", "VitA+MeHg"),
  treatment_colors = c(
    CN = "#98F5FF",
    VitA = "#B3EE3A",
    MeHg = "#CD2626",
    `VitA+MeHg` = "#4D4D4D"
  ),
  expected_sample_counts = c(Female = 18L),
  expected_treatment_counts = list(
    Female = c(CN = 5L, VitA = 4L, MeHg = 5L, `VitA+MeHg` = 4L)
  ),
  composite_filename = "Figure_S2_FinalFigure.png",
  composite_dimensions = c(2082, 629),
  expected_abundance_md5 = "92beb6c9b419284b9dc8cd3574cd7a86",
  expected_manifest_md5 = "17734fb98113e31be5d2b8834e832b24",
  expected_omnibus_statistics = c(
    12.045260347129489,
    10.792697290930505,
    8.611007462686564
  ),
  expected_omnibus_p = c(
    0.007229714926238264,
    0.012901315864119795,
    0.034935804811206110
  ),
  expected_significant_counts = c(3L, 2L, 1L),
  expected_significant_signatures = c(
    "CN|MeHg,VitA|MeHg,MeHg|VitA+MeHg",
    "CN|VitA+MeHg,MeHg|VitA+MeHg",
    "MeHg|VitA+MeHg"
  ),
  expected_significant_adjusted_p = list(
    c(
      0.021116679185680144,
      0.013774304703013730,
      0.025317525090288416
    ),
    c(0.023790436661229153, 0.012728977803228965),
    0.049157091296109755
  )
)

generate_figure_s2_abc <- function(
    output_dir = file.path(
      find_supplementary_abundance_project_root(),
      "supplementary_figures",
      "Figure_S2",
      "output"
    )) {
  generate_supplementary_abundance_figure(
    figure_s2_specification,
    output_dir
  )
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s2_abc()
  message(
    "Generated Figure S2 panels A-C in: ",
    normalizePath(dirname(result$files[[1L]]), winslash = "/")
  )
}
