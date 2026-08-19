# Reproducible Figure S4 panels A-AN.
#
# Figure S4 contains the 29 female and 11 male F1-PND365 taxon-abundance
# panels used in the supplied Photoshop composite. This script generates the
# individual panels only.

find_figure_s4_common_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      file.path(path, "figure_s3_s5_abundance_common.R"),
      file.path(
        path,
        "supplementary_figures",
        "figure_s3_s5_abundance_common.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate figure_s3_s5_abundance_common.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s4_common_script())

figure_s4_labels <- c(
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
  "g__.Eubacterium..ruminantium.group",
  "g__Alistipes",
  "g__Allobaculum",
  "g__Amylolactobacillus",
  "g__ASF356",
  "g__Bifidobacterium",
  "g__CAG.352",
  "g__Christensenellaceae.R.7.group",
  "g__Dubosiella",
  "g__HT002",
  "g__Lachnospiraceae.NK4A136.group",
  "g__Limosilactobacillus",
  "g__Phascolarctobacterium",
  "g__Prevotellaceae.Ga6A1.group",
  "c__Alphaproteobacteria",
  "o__Clostridiales",
  "o__Rickettsiales",
  "f__Clostridiaceae",
  "f__Mitochondria",
  "g__Alistipes",
  "g__CAG.352",
  "g__Clostridium.sensu.stricto.2",
  "g__Lactovum",
  "g__Prevotella",
  "g__Streptococcus"
)

figure_s4_specification <- list(
  figure_id = "S4",
  batch_id = unname(supplementary_abundance_batch_ids[["batch2"]]),
  generation = "F1",
  panels = data.frame(
    panel = supplementary_abundance_panel_codes(length(figure_s4_labels)),
    sex = c(rep("Female", 29L), rep("Male", 11L)),
    label = figure_s4_labels,
    stringsAsFactors = FALSE
  ),
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
  composite_filename = "Figure_S4_FinalFigure.png",
  composite_dimensions = c(4668, 5929),
  expected_abundance_md5 = "f3f65293ad4e35e647e6d93d393c5afb",
  expected_manifest_md5 = "a703799a2d60b1c3d935901ad76d8d6c",
  expected_omnibus_statistics = c(
    9.56641604010025, 9.56641604010025, 9.46153846153846,
    8.72826086956522, 8.72826086956522, 9.56641604010025,
    8.74358974358974, 8.43589743589744, 8.12820512820513,
    8.72826086956522, 9.39181286549707, 9.56641604010025,
    8.74358974358974, 8.43589743589744, 8.23076923076923,
    8.74154589371981, 9.39181286549707, 9.05128205128205,
    9.31188118811882, 8.19420289855073, 9.56641604010025,
    8.07835325365206, 8.74358974358974, 10.7349397590362,
    8.07692307692308, 7.82051282051282, 8.19420289855073,
    8.72826086956522, 9.18045112781955, 8.33333333333334,
    9.31188118811882, 8.33333333333334, 9.31188118811882,
    8.33333333333334, 8.30314009661836, 9.31188118811882,
    9.31188118811882, 8.29537953795379, 8.12820512820513,
    8.43589743589744
  ),
  expected_omnibus_p = c(
    0.022635201644082, 0.022635201644082, 0.0237440664767199,
    0.0331307275220085, 0.0331307275220085, 0.022635201644082,
    0.0329015914281892, 0.0378118050742007, 0.0434350841537281,
    0.0331307275220085, 0.0245105839366782, 0.022635201644082,
    0.0329015914281892, 0.0378118050742007, 0.0414755598461637,
    0.032932053125349, 0.0245105839366782, 0.028616776874436,
    0.0254191830349234, 0.0421640767150978, 0.022635201644082,
    0.0444198331305124, 0.0329015914281892, 0.0132489965279209,
    0.044448402794751, 0.0498705297687475, 0.0421640767150978,
    0.0331307275220085, 0.0269854524770431, 0.0396023552075641,
    0.0254191830349234, 0.0396023552075641, 0.0254191830349234,
    0.0396023552075641, 0.0401450508018971, 0.0254191830349234,
    0.0254191830349234, 0.0402857078711178, 0.0434350841537281,
    0.0378118050742007
  ),
  expected_significant_counts = c(
    1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 2L, 1L, 0L, 0L, 0L,
    0L, 2L, 0L, 2L, 1L, 1L, 1L, 0L, 3L, 1L, 1L, 2L, 1L, 0L, 0L,
    2L, 0L, 2L, 0L, 1L, 2L, 2L, 2L, 1L, 0L
  ),
  expected_significant_signatures = c(
    "CN|VitA+MeHg", "CN|VitA+MeHg", "VitA|VitA+MeHg", "CN|MeHg",
    "CN|MeHg", "CN|VitA+MeHg", "", "", "", "CN|MeHg",
    "CN|MeHg,VitA|MeHg", "CN|VitA+MeHg", "", "", "", "",
    "CN|MeHg,VitA|MeHg", "", "VitA|MeHg,VitA|VitA+MeHg",
    "MeHg|VitA+MeHg", "CN|VitA+MeHg", "VitA|MeHg", "",
    "CN|VitA,CN|MeHg,CN|VitA+MeHg", "VitA|MeHg", "MeHg|VitA+MeHg",
    "VitA|MeHg,VitA|VitA+MeHg", "CN|MeHg", "", "",
    "CN|VitA,CN|MeHg", "", "CN|VitA,CN|MeHg", "",
    "VitA|VitA+MeHg", "CN|MeHg,CN|VitA+MeHg", "CN|VitA,CN|MeHg",
    "CN|VitA,CN|VitA+MeHg", "CN|MeHg", ""
  )
)

generate_figure_s4_a_an <- function(
    output_dir = file.path(
      find_supplementary_abundance_project_root(),
      "supplementary_figures",
      "Figure_S4",
      "output"
    )) {
  generate_supplementary_abundance_figure(
    figure_s4_specification,
    output_dir
  )
}

if (sys.nframe() == 0L) {
  result <- generate_figure_s4_a_an()
  message(
    "Generated Figure S4 panels A-AN in: ",
    normalizePath(dirname(result$files[[1L]]), winslash = "/")
  )
}
