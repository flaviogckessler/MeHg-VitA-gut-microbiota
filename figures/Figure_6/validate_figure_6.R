# Validate Figure 6 sources, cohort selection, provenance, scripts, and PNGs.
# Set FIGURE6_FULL_VALIDATION=true to rerun the three seeded ALDEx2 generators;
# the default structural validation is fast because each full panel rerun can
# take several minutes on Windows.

# Load Packages ----
suppressPackageStartupMessages(library(utils))

find_figure_6_dir <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    candidates <- c(path, file.path(path, "figures", "Figure_6"))
    for (candidate in candidates) {
      required <- file.path(
        candidate,
        c(
          "figure_6_common.R",
          "generate_figure_6_a.R",
          "generate_figure_6_b.R",
          "generate_figure_6_c.R"
        )
      )
      if (all(file.exists(required))) return(candidate)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate the Figure 6 workflow.", call. = FALSE)
    }
    path <- parent
  }
}

check_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_png_dimensions <- function(filename) {
  connection <- file(filename, open = "rb")
  on.exit(close(connection), add = TRUE)
  header <- readBin(connection, what = "raw", n = 24L)
  signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(header) < 24L || !identical(header[1:8], signature)) {
    stop("Invalid PNG signature: ", filename, call. = FALSE)
  }
  unsigned_32 <- function(bytes) {
    sum(as.numeric(as.integer(bytes)) * 256^(3:0))
  }
  c(width = unsigned_32(header[17:20]), height = unsigned_32(header[21:24]))
}

validate_figure_6 <- function(
    rerun = identical(tolower(Sys.getenv("FIGURE6_FULL_VALIDATION")), "true")) {
  figure_dir <- find_figure_6_dir()
  project_root <- normalizePath(
    file.path(figure_dir, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  scripts <- file.path(
    figure_dir,
    c(
      "figure_6_common.R",
      "generate_figure_6_a.R",
      "generate_figure_6_b.R",
      "generate_figure_6_c.R"
    )
  )
  invisible(lapply(scripts, parse))

  batch_ids <- c(
    "batch01_iontorrent_dams-f0_f1-pnd21",
    "batch02_illumina-miseq_f1-f2_pnd365"
  )
  ko_files <- file.path(
    project_root,
    "results", "functional_prediction", batch_ids,
    "kegg_ortholog_abundance.tsv"
  )
  manifest_file <- file.path(
    project_root, "data", "metadata", "sample_manifest_all_batches.csv"
  )
  canonical_inputs <- c(ko_files, manifest_file)
  check_true(
    all(file.exists(canonical_inputs)) && all(file.info(canonical_inputs)$size > 0),
    "One or more canonical Figure 6 inputs are missing or empty."
  )

  first_ids <- vapply(ko_files, function(filename) {
    table <- utils::read.delim(
      filename,
      nrows = 1L,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
    table$ID[[1L]]
  }, character(1))
  check_true(
    all(grepl("^ko:K[0-9]{5}$", first_ids)),
    "Canonical KO inputs do not retain the expected PICRUSt2 ko:K##### IDs."
  )

  manifest <- utils::read.csv(
    manifest_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM"
  )
  included <- tolower(as.character(manifest$include_primary_analysis)) %in%
    c("true", "t", "1", "yes")
  strata <- list(
    A = manifest[
      manifest$batch_id == batch_ids[[1L]] & manifest$generation == "F1" & included,
      ,
      drop = FALSE
    ],
    B = manifest[
      manifest$batch_id == batch_ids[[2L]] & manifest$generation == "F1" & included,
      ,
      drop = FALSE
    ],
    C = manifest[
      manifest$batch_id == batch_ids[[2L]] & manifest$generation == "F2" & included,
      ,
      drop = FALSE
    ]
  )
  check_true(
    identical(vapply(strata, nrow, integer(1)), c(A = 34L, B = 24L, C = 24L)),
    "Figure 6 primary-analysis sample counts are incorrect."
  )
  check_true(
    !any(c("10F", "22F", "t3") %in% unlist(lapply(strata, `[[`, "sample_id"))),
    "An excluded or undocumented sample entered a Figure 6 stratum."
  )
  for (panel in names(strata)) {
    counts <- table(factor(
      strata[[panel]]$treatment,
      levels = c("CN", "VitA", "MeHg", "VitA+MeHg")
    ))
    check_true(
      all(counts > 0L),
      paste("Panel", panel, "does not contain all four treatment groups.")
    )
  }

  legacy_file <- file.path(
    project_root, "scripts", "legacy_scripts", "Picrust_mestrado_Alberto.R"
  )
  old_legacy_location <- file.path(project_root, "scripts", "Picrust_mestrado_Alberto.R")
  check_true(
    file.exists(legacy_file) && !file.exists(old_legacy_location),
    "The Figure 6 legacy reference is not isolated under scripts/legacy_scripts/."
  )

  output_files <- file.path(
    figure_dir,
    "output",
    c(
      "Figure6A_generated.png",
      "Figure6B_generated.png",
      "Figure6C_generated.png"
    )
  )
  expected_dimensions <- list(
    c(width = 2048, height = 760),
    c(width = 2048, height = 265),
    c(width = 2048, height = 265)
  )
  check_true(
    all(file.exists(output_files)) && all(file.info(output_files)$size > 5000),
    "One or more generated Figure 6 PNGs are missing or empty."
  )
  for (index in seq_along(output_files)) {
    check_true(
      identical(read_png_dimensions(output_files[[index]]), expected_dimensions[[index]]),
      paste("Generated panel", LETTERS[[index]], "has incorrect PNG dimensions.")
    )
  }
  reference_file <- file.path(figure_dir, "output", "Figure_6_FinalFigure.png")
  check_true(
    file.exists(reference_file) &&
      identical(read_png_dimensions(reference_file), c(width = 2198, height = 1424)),
    "The supplied Figure 6 reference PNG is missing or has changed dimensions."
  )

  if (isTRUE(rerun)) {
    panel_letters <- c("a", "b", "c")
    environments <- lapply(panel_letters, function(panel) {
      environment <- new.env(parent = globalenv())
      sys.source(
        file.path(figure_dir, paste0("generate_figure_6_", panel, ".R")),
        envir = environment
      )
      environment
    })
    validation_output <- file.path(tempdir(), "figure_6_full_validation")
    results <- Map(function(environment, panel) {
      generator <- get(paste0("generate_figure_6_", panel), envir = environment)
      generator(output_dir = validation_output)
    }, environments, panel_letters)
    expected_features <- list(
      c(
        "ko03050", "ko03450", "ko00310", "ko00190", "ko00740",
        "ko00760", "ko00430", "ko00480", "ko00623", "ko00930"
      ),
      "ko03050",
      "ko00311"
    )
    for (index in seq_along(results)) {
      check_true(
        setequal(results[[index]]$selected_daa$feature, expected_features[[index]]),
        paste("Full validation selected incorrect pathways for panel", LETTERS[[index]], ".")
      )
    }
  }

  message(
    "Figure 6 passed source, cohort, provenance, script, reference-image, ",
    "and PNG-output validation",
    if (isTRUE(rerun)) " including full seeded ALDEx2 reruns." else "."
  )
  invisible(TRUE)
}

if (sys.nframe() == 0L) validate_figure_6()
