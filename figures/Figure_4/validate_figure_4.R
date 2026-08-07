# Validate the four reproducible Figure 4 panel workflows.
# This follows the Figure 3 validation pattern while also checking the preserved
# significant-taxon checkpoints under results/legacy_summary_tables/figure_4/.

# - Locate and load each panel script in an isolated environment
find_figure_4_dir <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(
      path,
      file.path(path, "figures", "Figure_4")
    )
    for (candidate in candidates) {
      required <- file.path(
        candidate,
        paste0("generate_figure_4_", letters[1:4], ".R")
      )
      if (all(file.exists(required))) return(candidate)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the Figure 4 generation scripts.", call. = FALSE)
    }
    path <- parent
  }
}

check_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

figure_dir <- find_figure_4_dir()
panel_names <- LETTERS[1:4]
script_paths <- setNames(
  file.path(figure_dir, paste0("generate_figure_4_", letters[1:4], ".R")),
  panel_names
)
panel_environments <- lapply(script_paths, function(script) {
  environment <- new.env(parent = globalenv())
  sys.source(script, envir = environment)
  environment
})

project_root <- panel_environments$A$find_project_root()

# - Confirm all runtime inputs are available inside New_repo
canonical_inputs <- c(
  batch1_abundance = file.path(
    project_root, "results", "taxon_relative_abundance",
    "batch01_iontorrent_dams-f0_f1-pnd21_taxon_relative_abundance.csv"
  ),
  batch2_abundance = file.path(
    project_root, "results", "taxon_relative_abundance",
    "batch02_illumina-miseq_f1-f2_pnd365_taxon_relative_abundance.csv"
  ),
  batch1_manifest = file.path(
    project_root, "data", "metadata",
    "batch01_iontorrent_dams-f0_f1-pnd21", "sample_manifest.csv"
  ),
  batch2_manifest = file.path(
    project_root, "data", "metadata",
    "batch02_illumina-miseq_f1-f2_pnd365", "sample_manifest.csv"
  )
)
check_true(
  all(file.exists(canonical_inputs)),
  "One or more canonical Figure 4 inputs are missing from New_repo."
)
check_true(
  all(file.info(canonical_inputs)$size > 0),
  "One or more canonical Figure 4 inputs are empty."
)

# - Prepare panel data
panel_data <- list(
  A = panel_environments$A$prepare_figure_4_a_data(project_root),
  B = panel_environments$B$prepare_figure_4_b_data(project_root),
  C = panel_environments$C$prepare_figure_4_c_data(project_root),
  D = panel_environments$D$prepare_figure_4_d_data(project_root)
)

expected_taxon_counts <- c(A = 3L, B = 17L, C = 40L, D = 15L)
observed_taxon_counts <- vapply(panel_data, nrow, integer(1))
check_true(
  identical(observed_taxon_counts, expected_taxon_counts),
  "One or more panels do not contain the expected number of significant taxa."
)

expected_stratum_counts <- list(
  A = c(Dams = 3L),
  B = c(Female = 10L, Male = 7L),
  C = c(Female = 29L, Male = 11L),
  D = c(Female = 8L, Male = 7L)
)
for (panel in panel_names) {
  observed <- table(factor(
    panel_data[[panel]]$stratum,
    levels = names(expected_stratum_counts[[panel]])
  ))
  check_true(
    identical(as.integer(observed), unname(expected_stratum_counts[[panel]])),
    paste("Panel", panel, "stratum row counts are incorrect.")
  )
}

expected_sample_counts <- list(
  A = c(Dams = 18L),
  B = c(Female = 18L, Male = 16L),
  C = c(Female = 12L, Male = 12L),
  D = c(Female = 12L, Male = 12L)
)
check_true(
  length(attr(panel_data$A, "sample_ids")) == expected_sample_counts$A[[1]],
  "Panel A must analyze 18 included dams."
)
for (panel in c("B", "C", "D")) {
  observed <- vapply(attr(panel_data[[panel]], "sample_ids"), length, integer(1))
  check_true(
    identical(observed, expected_sample_counts[[panel]]),
    paste("Panel", panel, "sample counts are incorrect.")
  )
}
check_true(
  !any(c("10F", "22F") %in% attr(panel_data$A, "sample_ids")),
  "Panel A contains batch-1 samples excluded for low sequencing quality."
)

# - Statistical and z-score validation
for (panel in panel_names) {
  data <- panel_data[[panel]]
  z_scores <- as.matrix(data[, paste0("z_", panel_environments$A$treatment_levels)])
  storage.mode(z_scores) <- "double"
  check_true(
    all(is.finite(data$p_value)) && all(data$p_value <= 0.05),
    paste("Panel", panel, "contains an invalid or nonsignificant p-value.")
  )
  check_true(
    all(abs(rowMeans(z_scores)) < 1e-12),
    paste("Panel", panel, "z-score rows do not have mean zero.")
  )
  check_true(
    all(abs(apply(z_scores, 1, stats::sd) - 1) < 1e-12),
    paste("Panel", panel, "z-score rows do not have unit sample SD.")
  )
  check_true(
    all(z_scores >= -1.5 - 1e-12 & z_scores <= 1.5 + 1e-12),
    paste("Panel", panel, "z-scores exceed the final legend limits.")
  )
}

# - Compare exact taxa, group means, and z-scores with legacy checkpoints
legacy_checkpoint_map <- list(
  A = c(Dams = "Signif_dams.csv"),
  B = c(
    Female = "Signif_offs_F1Female.csv",
    Male = "Signif_offs_F1Male.csv"
  ),
  C = c(
    Female = "Signif_1year_F1Female.csv",
    Male = "Signif_1year_F1Male.csv"
  ),
  D = c(
    Female = "Signif_1year_F2Female.csv",
    Male = "Signif_1year_F2Male.csv"
  )
)
legacy_dir <- file.path(
  project_root, "results", "legacy_summary_tables", "figure_4"
)
treatments <- panel_environments$A$treatment_levels

for (panel in panel_names) {
  for (stratum in names(legacy_checkpoint_map[[panel]])) {
    checkpoint_file <- file.path(
      legacy_dir, legacy_checkpoint_map[[panel]][[stratum]]
    )
    checkpoint <- read.csv(
      checkpoint_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    observed <- panel_data[[panel]][
      panel_data[[panel]]$stratum == stratum,
      ,
      drop = FALSE
    ]
    check_true(
      identical(observed$taxon, checkpoint$tax),
      paste("Panel", panel, stratum, "taxon order differs from the legacy checkpoint.")
    )

    observed_means <- as.matrix(
      observed[, paste0("mean_", treatments), drop = FALSE]
    )
    checkpoint_means <- as.matrix(checkpoint[, treatments, drop = FALSE])
    storage.mode(observed_means) <- "double"
    storage.mode(checkpoint_means) <- "double"
    check_true(
      max(abs(observed_means - checkpoint_means)) < 1e-12,
      paste("Panel", panel, stratum, "group means differ from the legacy checkpoint.")
    )

    checkpoint_z <- t(apply(checkpoint_means, 1, function(x) {
      (x - mean(x)) / stats::sd(x)
    }))
    observed_z <- as.matrix(observed[, paste0("z_", treatments), drop = FALSE])
    storage.mode(observed_z) <- "double"
    check_true(
      max(abs(observed_z - checkpoint_z)) < 1e-12,
      paste("Panel", panel, stratum, "z-scores differ from the legacy checkpoint.")
    )
  }
}

# - PNG helper and reference-image checks
read_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)
  header <- readBin(connection, what = "raw", n = 24L)
  signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(header) < 24L || !identical(header[1:8], signature)) {
    stop(paste("Invalid PNG signature:", filename), call. = FALSE)
  }
  unsigned_32 <- function(bytes) {
    sum(as.numeric(as.integer(bytes)) * 256^(3:0))
  }
  c(width = unsigned_32(header[17:20]), height = unsigned_32(header[21:24]))
}

reference_file <- file.path(figure_dir, "output", "Figure_4_FinalFigure.png")
check_true(file.exists(reference_file), "The final Figure 4 reference PNG is missing.")
check_true(
  identical(unname(read_png_dimensions(reference_file)), c(2520, 1514)),
  "The final Figure 4 reference dimensions have changed."
)

# - Output validation in an isolated temporary directory
validation_dir <- file.path(tempdir(), "figure_4_validation")
generation_results <- list(
  A = panel_environments$A$generate_figure_4_a(validation_dir),
  B = panel_environments$B$generate_figure_4_b(validation_dir),
  C = panel_environments$C$generate_figure_4_c(validation_dir),
  D = panel_environments$D$generate_figure_4_d(validation_dir)
)
generated_files <- vapply(generation_results, `[[`, character(1), "file")
expected_names <- paste0("Figure4", panel_names, "_generated.png")
expected_dimensions <- list(
  A = c(width = 1000, height = 360),
  B = c(width = 1150, height = 825),
  C = c(width = 1200, height = 1532),
  D = c(width = 1250, height = 810)
)

check_true(all(file.exists(generated_files)), "One or more Figure 4 PNGs were not generated.")
check_true(
  identical(basename(generated_files), expected_names),
  "Generated Figure 4 panel filenames do not match the expected names."
)
check_true(
  all(file.info(generated_files)$size > 1000),
  "One or more Figure 4 PNGs are empty."
)
for (panel in panel_names) {
  check_true(
    identical(read_png_dimensions(generated_files[[panel]]), expected_dimensions[[panel]]),
    paste("Panel", panel, "PNG dimensions are incorrect.")
  )
}

message(
  "Figure 4 panels A-D passed source, statistics, legacy-checkpoint, ",
  "reference-image, and PNG-output validation."
)
