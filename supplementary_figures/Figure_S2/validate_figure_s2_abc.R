# Validation for the Figure S2 panels A-C workflow.

find_figure_s2_script <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    candidates <- c(
      file.path(path, "generate_figure_s2_abc.R"),
      file.path(
        path,
        "supplementary_figures",
        "Figure_S2",
        "generate_figure_s2_abc.R"
      )
    )
    existing <- candidates[file.exists(candidates)]
    if (length(existing) > 0L) return(existing[[1L]])

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate generate_figure_s2_abc.R.", call. = FALSE)
    }
    path <- parent
  }
}

source(find_figure_s2_script())

check_figure_s2_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

check_figure_s2_numeric <- function(observed, expected, tolerance, message) {
  check_figure_s2_true(
    length(observed) == length(expected) &&
      all(is.finite(observed)) &&
      isTRUE(all.equal(
        as.numeric(observed),
        as.numeric(expected),
        tolerance = tolerance,
        check.attributes = FALSE
      )),
    message
  )
}

read_figure_s2_png_dimensions <- function(filename) {
  header <- readBin(filename, what = "raw", n = 24L)
  png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (length(header) != 24L || !identical(header[1:8], png_signature)) {
    stop("Invalid PNG file: ", filename, call. = FALSE)
  }
  uint32 <- function(bytes) {
    sum(as.integer(bytes) * 256^(3:0))
  }
  c(
    width = uint32(header[17:20]),
    height = uint32(header[21:24])
  )
}

# Input and cohort checks
project_root <- find_figure_s2_project_root()
sources <- read_figure_s2_sources(project_root)
panel_data <- prepare_figure_s2_panel_data(sources)

check_figure_s2_true(
  identical(names(panel_data), c("A", "B", "C")),
  "Figure S2 panel IDs must be A-C."
)
check_figure_s2_true(
  all(vapply(panel_data, nrow, integer(1)) == 18L),
  "Every Figure S2 panel must contain 18 F0 dams."
)
check_figure_s2_true(
  all(vapply(panel_data, function(data) {
    identical(
      as.integer(table(data$treatment)),
      c(5L, 4L, 5L, 4L)
    ) && identical(names(table(data$treatment)), figure_s2_treatment_levels)
  }, logical(1))),
  "Figure S2 treatment counts or ordering changed."
)
check_figure_s2_true(
  all(vapply(panel_data, function(data) {
    length(unique(data$sample_id)) == 18L &&
      all(is.finite(data$relative_abundance)) &&
      all(data$relative_abundance >= 0 & data$relative_abundance <= 1)
  }, logical(1))),
  "Figure S2 contains duplicated samples or invalid abundances."
)
check_figure_s2_true(
  all(vapply(figure_s2_panel_configurations, function(configuration) {
    identical(make.names(configuration$taxon), configuration$label)
  }, logical(1))),
  "Figure S2 display labels no longer reproduce the legacy make.names labels."
)

# Nonparametric statistical checks
statistics <- lapply(panel_data, calculate_figure_s2_statistics)
expected_omnibus_statistics <- c(
  A = 12.045260347129489,
  B = 10.792697290930505,
  C = 8.611007462686564
)
expected_omnibus_p <- c(
  A = 0.007229714926238264,
  B = 0.012901315864119795,
  C = 0.034935804811206110
)
observed_omnibus_statistics <- vapply(
  statistics,
  function(result) result$omnibus$statistic[[1L]],
  numeric(1)
)
observed_omnibus_p <- vapply(
  statistics,
  function(result) result$omnibus$p[[1L]],
  numeric(1)
)

check_figure_s2_numeric(
  observed_omnibus_statistics,
  expected_omnibus_statistics,
  tolerance = 1e-12,
  "Figure S2 Kruskal-Wallis statistics changed."
)
check_figure_s2_numeric(
  observed_omnibus_p,
  expected_omnibus_p,
  tolerance = 1e-12,
  "Figure S2 Kruskal-Wallis p-values changed."
)
check_figure_s2_true(
  all(vapply(statistics, function(result) {
    identical(result$omnibus$method, "Kruskal-Wallis") &&
      result$omnibus$n[[1L]] == 18L &&
      result$omnibus$df[[1L]] == 3L &&
      nrow(result$pairwise) == 6L
  }, logical(1))),
  "Figure S2 must use Kruskal-Wallis tests and all six Dunn comparisons."
)

expected_significant_pairs <- list(
  A = c("CN|MeHg", "VitA|MeHg", "MeHg|VitA+MeHg"),
  B = c("CN|VitA+MeHg", "MeHg|VitA+MeHg"),
  C = "MeHg|VitA+MeHg"
)
expected_adjusted_p <- list(
  A = c(
    0.021116679185680144,
    0.013774304703013730,
    0.025317525090288416
  ),
  B = c(0.023790436661229153, 0.012728977803228965),
  C = 0.049157091296109755
)
for (panel in names(statistics)) {
  significant <- statistics[[panel]]$significant
  observed_pairs <- paste(significant$group1, significant$group2, sep = "|")
  check_figure_s2_true(
    identical(observed_pairs, expected_significant_pairs[[panel]]),
    paste("Figure S2 panel", panel, "significant Dunn pairs changed.")
  )
  check_figure_s2_numeric(
    significant$p.adj,
    expected_adjusted_p[[panel]],
    tolerance = 1e-12,
    paste("Figure S2 panel", panel, "FDR-adjusted Dunn p-values changed.")
  )
  check_figure_s2_true(
    all(significant$p.adj.signif == "*") &&
      all(diff(significant$y.position) > 0),
    paste("Figure S2 panel", panel, "annotation ordering changed.")
  )
}

# Rendering checks
validation_dir <- file.path(tempdir(), "figure_s2_abc_validation")
rendered <- generate_figure_s2_abc(validation_dir)
expected_names <- paste0("FigureS2", LETTERS[1:3], "_generated.png")
check_figure_s2_true(
  identical(basename(unname(rendered$files)), expected_names),
  "Figure S2 generated filenames changed."
)
check_figure_s2_true(
  all(file.exists(rendered$files)),
  "One or more Figure S2 validation PNGs were not generated."
)
check_figure_s2_true(
  all(vapply(rendered$files, function(filename) {
    identical(
      unname(read_figure_s2_png_dimensions(filename)),
      c(600, 500)
    )
  }, logical(1))),
  "Figure S2 validation PNGs must be 600 x 500 pixels."
)
check_figure_s2_true(
  all(file.info(rendered$files)$size > 5000 &
        file.info(rendered$files)$size < 100000),
  "A Figure S2 validation PNG is empty or unexpectedly large."
)
check_figure_s2_true(
  all(vapply(rendered$plots, inherits, logical(1), what = "ggplot")),
  "One or more Figure S2 panel objects are not ggplot objects."
)

committed_files <- file.path(
  project_root,
  "supplementary_figures",
  "Figure_S2",
  "output",
  expected_names
)
check_figure_s2_true(
  all(file.exists(committed_files)) &&
    all(vapply(committed_files, function(filename) {
      identical(
        unname(read_figure_s2_png_dimensions(filename)),
        c(600, 500)
      )
    }, logical(1))),
  "Committed Figure S2 generated panels are missing or have invalid dimensions."
)

reference_composite <- file.path(
  project_root,
  "supplementary_figures",
  "Figure_S2",
  "output",
  "Figure_S2_FinalFigure.png"
)
check_figure_s2_true(
  file.exists(reference_composite) &&
    identical(
      unname(read_figure_s2_png_dimensions(reference_composite)),
      c(2082, 629)
    ),
  "The supplied Figure S2 reference composite is missing or has changed dimensions."
)

message(
  "Figure S2 panels A-C passed source, cohort, nonparametric-statistics, ",
  "annotation, and PNG-rendering checks."
)
