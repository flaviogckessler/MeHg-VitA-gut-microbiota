# Generate one accessible canonical SVG from workflow_spec.json, rasterize that
# SVG to PNG, validate staged artifacts, and only then promote canonical files.
# Scientific wording, mappings, decisions, styles, colors, and layout are JSON-owned.

locate_workflow_directory <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) > 0L) {
    script_path <- sub("^--file=", "", file_argument[[1L]])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }
  candidates <- c(file.path(getwd(), "figures", "Figure_1", "Workflow"), getwd())
  matches <- candidates[file.exists(file.path(candidates, "workflow_spec_lib.R"))]
  if (length(matches) == 0L) {
    stop("Could not locate figures/Figure_1/Workflow/workflow_spec_lib.R.", call. = FALSE)
  }
  normalizePath(matches[[1L]], winslash = "/", mustWork = TRUE)
}

workflow_directory <- locate_workflow_directory()
source(file.path(workflow_directory, "workflow_spec_lib.R"), local = TRUE)

xml_escape <- function(value) {
  value <- gsub("&", "&amp;", as.character(value), fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  gsub("'", "&apos;", value, fixed = TRUE)
}

palette_color <- function(spec, role) {
  value <- spec$render$theme$palette[[role]]
  if (is.null(value)) stop(sprintf("Unknown palette role: %s", role), call. = FALSE)
  value
}

node_index <- function(spec) workflow_index_by_id(spec$nodes)
card_index <- function(spec) workflow_index_by_id(spec$render$cards)
lane_index <- function(spec) workflow_index_by_id(spec$lanes)
style_index <- function(spec) workflow_index_by_id(spec$render$execution_styles)
lane_layout_index <- function(spec) {
  stats::setNames(spec$render$lane_layout, vapply(
    spec$render$lane_layout, function(item) item$lane_id, character(1)
  ))
}

card_title <- function(spec, card) {
  nodes <- node_index(spec)
  labels <- vapply(card$title$node_ids, function(id) nodes[[id]]$display_label, character(1))
  paste(labels, collapse = card$title$separator)
}

card_items <- function(spec, card) {
  nodes <- node_index(spec)
  items <- list()
  for (source in card$bullet_sources) {
    values <- unlist(nodes[[source$node_id]][[source$field]], use.names = FALSE)
    selected <- values[unlist(source$indexes, use.names = FALSE)]
    terms <- unique(c(
      unlist(nodes[[source$node_id]]$italic_terms, use.names = FALSE),
      unlist(source$italic_terms, use.names = FALSE)
    ))
    for (value in selected) {
      items[[length(items) + 1L]] <- list(text = value, italic_terms = terms)
    }
  }
  items
}

wrap_for_width <- function(text, available_width, font_size, bold = FALSE) {
  width_factor <- if (isTRUE(bold)) 0.54 else 0.49
  characters <- max(8L, floor(available_width / (font_size * width_factor)))
  wrapped <- strwrap(text, width = characters, simplify = TRUE)
  if (length(wrapped) == 0L) "" else wrapped
}

styled_svg_text <- function(text, italic_terms = character()) {
  if (length(italic_terms) == 0L) return(xml_escape(text))
  remaining <- text
  output <- character()
  while (nzchar(remaining)) {
    positions <- vapply(italic_terms, function(term) {
      match <- regexpr(term, remaining, fixed = TRUE)[[1L]]
      if (match < 0L) Inf else match
    }, numeric(1))
    if (all(is.infinite(positions))) {
      output <- c(output, xml_escape(remaining))
      break
    }
    selected <- which.min(positions)
    position <- positions[[selected]]
    term <- italic_terms[[selected]]
    if (position > 1L) {
      output <- c(output, xml_escape(substr(remaining, 1L, position - 1L)))
    }
    output <- c(
      output,
      sprintf('<tspan font-style="italic">%s</tspan>', xml_escape(term))
    )
    remaining <- substr(remaining, position + nchar(term), nchar(remaining))
  }
  paste0(output, collapse = "")
}

text_element <- function(class_name, role, x, y, font_size, lines,
                         line_height, weight = 400, italic_terms = character()) {
  tspans <- vapply(seq_along(lines), function(i) {
    baseline <- y + (i - 1L) * line_height
    sprintf(
      '<tspan x="%.2f" y="%.2f">%s</tspan>',
      x, baseline, styled_svg_text(lines[[i]], italic_terms)
    )
  }, character(1))
  sprintf(
    '<text class="%s" data-role="%s" x="%.2f" y="%.2f" font-size="%.2f" font-weight="%d">%s</text>',
    class_name, role, x, y, font_size, weight, paste0(tspans, collapse = "")
  )
}

lane_svg <- function(spec, layout) {
  lane <- lane_index(spec)[[layout$lane_id]]
  c(
    sprintf('<g id="lane-%s" class="workflow-lane">', xml_escape(lane$id)),
    sprintf(
      '<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="18" fill="%s" stroke="%s" stroke-width="1.5"/>',
      layout$x, layout$y, layout$width, layout$height,
      lane$background_color, lane$figure_color
    ),
    text_element(
      "lane-label", "lane-label",
      layout$x + 22, layout$y + 39, 30,
      lane$display_label, 30, 700
    ),
    sprintf(
      '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="4"/>',
      layout$x + 22, layout$y + 52, layout$x + 245, layout$y + 52,
      lane$figure_color
    ),
    '</g>'
  )
}

card_svg <- function(spec, card) {
  nodes <- node_index(spec)
  node <- nodes[[card$member_node_ids[[1L]]]]
  style <- style_index(spec)[[node$execution$visual_class]]
  accent <- if (identical(card$accent_role, "lane")) {
    lane_index(spec)[[card$lane_id]]$figure_color
  } else {
    palette_color(spec, card$accent_role)
  }
  dash <- if (is.null(style$stroke_dasharray)) "" else {
    sprintf(' stroke-dasharray="%s"', xml_escape(style$stroke_dasharray))
  }
  padding <- 16
  title_size <- 30
  body_size <- 30
  title_lines <- wrap_for_width(
    card_title(spec, card), card$width - 2 * padding, title_size, bold = TRUE
  )
  title_y <- card$y + 53
  title_line_height <- 29
  body_y <- title_y + (length(title_lines) - 1L) * title_line_height + 34

  body_lines <- character()
  body_terms <- character()
  for (item in card_items(spec, card)) {
    wrapped <- wrap_for_width(
      item$text, card$width - 2 * padding - 12, body_size, bold = FALSE
    )
    for (i in seq_along(wrapped)) {
      body_lines <- c(
        body_lines,
        paste0(if (i == 1L) "\u2022 " else "  ", wrapped[[i]])
      )
      body_terms <- unique(c(body_terms, item$italic_terms))
    }
  }

  node_title <- paste0(
    node$display_label, ". Execution: ", style$display_label,
    ". Evidence: ", node$evidence$status, " / ", node$evidence$authority, "."
  )
  c(
    sprintf(
      '<g id="node-%s" class="workflow-node execution-%s" data-card-id="%s" data-lane-id="%s">',
      xml_escape(node$id), xml_escape(node$execution$visual_class),
      xml_escape(card$id), xml_escape(card$lane_id)
    ),
    sprintf('<title>%s</title>', xml_escape(node_title)),
    sprintf(
      '<rect class="card-frame" x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="14" fill="%s" stroke="%s" stroke-width="%.2f"%s/>',
      card$x, card$y, card$width, card$height,
      palette_color(spec, card$fill_role), accent, style$stroke_width, dash
    ),
    sprintf(
      '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="6" stroke-linecap="round"/>',
      card$x + 12, card$y + 14, card$x + 12, card$y + card$height - 14, accent
    ),
    text_element(
      "card-status", "card-status",
      card$x + padding + 8, card$y + 24, 22,
      style$card_label, 22, 700
    ),
    text_element(
      "card-title", "card-title",
      card$x + padding + 8, title_y, title_size,
      title_lines, title_line_height, 700
    ),
    text_element(
      "card-body", "card-body",
      card$x + padding + 8, body_y, body_size,
      body_lines, 28, 400, body_terms
    ),
    '</g>'
  )
}

connector_segment_path <- function(segment) {
  points <- segment$points
  commands <- sprintf("M %.2f %.2f", points[[1L]]$x, points[[1L]]$y)
  for (i in 2:length(points)) {
    previous <- points[[i - 1L]]
    current <- points[[i]]
    if (isTRUE(all.equal(as.numeric(previous$y), as.numeric(current$y)))) {
      commands <- c(commands, sprintf("H %.2f", current$x))
    } else if (isTRUE(all.equal(as.numeric(previous$x), as.numeric(current$x)))) {
      commands <- c(commands, sprintf("V %.2f", current$y))
    } else {
      stop(sprintf("Connector segment is not orthogonal: %s", segment$id), call. = FALSE)
    }
  }
  paste(commands, collapse = " ")
}

connector_group_svg <- function(spec, group) {
  output <- sprintf(
    '<g id="connector-group-%s" class="connector-group" data-semantic-edge-ids="%s">',
    xml_escape(group$id),
    xml_escape(paste(unlist(group$semantic_edge_ids, use.names = FALSE), collapse = " "))
  )
  for (segment in group$segments) {
    marker <- if (isTRUE(segment$marker_end)) ' marker-end="url(#arrowhead)"' else ""
    output <- c(output, sprintf(
      '<path id="connector-%s-%s" class="workflow-connector" d="%s" fill="none" stroke="%s" stroke-width="%.2f" stroke-linecap="round" stroke-linejoin="round"%s/>',
      xml_escape(group$id), xml_escape(segment$id), connector_segment_path(segment),
      palette_color(spec, spec$render$connector_policy$color_role),
      spec$render$connector_policy$line_width, marker
    ))
  }
  c(output, '</g>')
}

legend_svg <- function(spec) {
  styles <- spec$render$execution_styles
  start_x <- 720
  item_width <- 255
  y <- 129
  output <- c('<g id="legend-execution" aria-label="Execution-style legend">')
  for (i in seq_along(styles)) {
    item <- styles[[i]]
    x <- start_x + (i - 1L) * item_width
    dash <- if (is.null(item$stroke_dasharray)) "" else {
      sprintf(' stroke-dasharray="%s"', xml_escape(item$stroke_dasharray))
    }
    output <- c(
      output,
      sprintf(
        '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%.2f"%s/>',
        x, y - 8, x + 42, y - 8, palette_color(spec, "ink"),
        item$stroke_width, dash
      ),
      text_element(
        "legend-label", "legend-label", x + 52, y, 22,
        item$card_label, 22, 700
      )
    )
  }
  c(output, '</g>')
}

svg_document <- function(spec) {
  canvas <- spec$render$canvas
  journal <- workflow_index_by_id(spec$render$profiles)[["journal_landscape"]]
  accessibility <- spec$render$accessibility
  font_stack <- paste(
    c(spec$render$font$family, unlist(spec$render$font$fallback_families)),
    collapse = ", "
  )
  metadata <- list(
    workflow_id = spec$workflow_id,
    schema_version = spec$schema_version,
    generator_version = spec$generator_contract$generator_version,
    profiles = spec$render$profiles
  )
  workflow_require_namespace("jsonlite")
  metadata_json <- jsonlite::toJSON(metadata, auto_unbox = TRUE, pretty = FALSE)

  output <- c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf(
      '<svg xmlns="http://www.w3.org/2000/svg" width="%smm" height="%smm" viewBox="0 0 %d %d" preserveAspectRatio="xMidYMid meet" role="img" aria-labelledby="%s %s" style="max-width:100%%;height:auto">',
      journal$width_mm, journal$height_mm,
      canvas$viewbox_width, canvas$viewbox_height,
      xml_escape(accessibility$title_id), xml_escape(accessibility$desc_id)
    ),
    sprintf('<title id="%s">%s</title>', xml_escape(accessibility$title_id), xml_escape(spec$title)),
    sprintf('<desc id="%s">%s</desc>', xml_escape(accessibility$desc_id), xml_escape(accessibility$alt_text)),
    sprintf('<metadata id="workflow-metadata">%s</metadata>', xml_escape(metadata_json)),
    '<defs>',
    sprintf(
      '<marker id="arrowhead" markerWidth="%d" markerHeight="%d" refX="8" refY="4.5" orient="auto" markerUnits="strokeWidth"><path d="M 0 0 L 9 4.5 L 0 9 z" fill="%s"/></marker>',
      spec$render$connector_policy$arrow_size,
      spec$render$connector_policy$arrow_size,
      palette_color(spec, spec$render$connector_policy$color_role)
    ),
    '</defs>',
    sprintf(
      '<style>text{font-family:%s;fill:%s}.card-status{fill:%s;letter-spacing:0.8px}.lane-label{fill:%s}.legend-label{fill:%s}</style>',
      xml_escape(font_stack), palette_color(spec, "text"),
      palette_color(spec, "muted"), palette_color(spec, "ink"),
      palette_color(spec, "text")
    ),
    sprintf(
      '<rect id="canvas-background" x="0" y="0" width="%d" height="%d" fill="%s"/>',
      canvas$viewbox_width, canvas$viewbox_height, palette_color(spec, "paper")
    ),
    text_element("figure-title", "figure-title", 42, 55, 38, spec$title, 38, 700),
    text_element("figure-subtitle", "figure-subtitle", 42, 101, 28, spec$render$subtitle, 28, 400),
    legend_svg(spec)
  )

  for (layout in spec$render$lane_layout) output <- c(output, lane_svg(spec, layout))
  output <- c(output, '<g id="workflow-connectors" aria-hidden="true">')
  for (group in spec$render$connector_groups) {
    output <- c(output, connector_group_svg(spec, group))
  }
  output <- c(output, '</g>', '<g id="workflow-nodes">')
  for (card in spec$render$cards) output <- c(output, card_svg(spec, card))
  c(output, '</g>', '</svg>')
}

caption_document <- function(spec) {
  lines <- c(
    "# Analysis workflow accessibility text",
    "",
    "## Alt text",
    "",
    spec$render$accessibility$alt_text,
    "",
    "## Caption",
    "",
    spec$render$caption$summary
  )
  if (isTRUE(spec$render$caption$include_exclusion_rules)) {
    lines <- c(lines, "", "### Sample exclusions", "")
    for (rule in spec$study_context$exclusion_rules) {
      lines <- c(lines, sprintf(
        "- `%s`: %s \u2014 %s (%s).",
        rule$scope,
        paste(unlist(rule$sample_ids), collapse = ", "),
        rule$reason,
        rule$rule
      ))
    }
  }
  if (isTRUE(spec$render$caption$include_manual_boundaries)) {
    lines <- c(lines, "", "### Scripted/manual boundaries", "")
    for (node in spec$nodes) {
      if (!is.null(node$execution$manual_boundary)) {
        lines <- c(lines, sprintf("- **%s:** %s", node$display_label, node$execution$manual_boundary))
      }
    }
  }
  lines <- c(lines, "", "### Figure notes", "")
  for (note in spec$render$caption$notes) lines <- c(lines, paste0("- ", note))
  paste(lines, collapse = "\n")
}

find_runtime <- function(environment_name, executable_name) {
  configured <- Sys.getenv(environment_name, unset = "")
  if (nzchar(configured)) return(configured)
  discovered <- Sys.which(executable_name)
  if (!nzchar(discovered)) {
    stop(sprintf(
      "%s is unavailable; set %s to its absolute path.",
      executable_name, environment_name
    ), call. = FALSE)
  }
  unname(discovered)
}

run_checked <- function(command, arguments, label) {
  quoted_arguments <- vapply(arguments, shQuote, character(1))
  output <- system2(command, quoted_arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(c(sprintf("%s failed:", label), output), collapse = "\n"), call. = FALSE)
  }
  if (length(output) > 0L) message(paste(output, collapse = "\n"))
  invisible(output)
}

write_utf8_lines <- function(lines, path) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(enc2utf8(lines), connection, sep = "\n", useBytes = TRUE)
}

write_staged_manifest <- function(spec, paths, staged) {
  workflow_require_namespace("jsonlite")
  identity <- workflow_freshness_identity(spec, paths)
  identity$generated_at_utc <- format(
    as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"
  )
  identity$profiles <- spec$render$profiles
  identity$font <- spec$render$font
  journal <- workflow_index_by_id(spec$render$profiles)[["journal_landscape"]]
  identity$outputs <- lapply(spec$generator_contract$outputs, function(output) {
    staged_path <- staged[[output$format]]
    record <- list(format = output$format, path = output$path, sha256 = workflow_sha256(staged_path))
    if (identical(output$format, "svg")) {
      record$width_mm <- journal$width_mm
      record$height_mm <- journal$height_mm
      record$viewbox_width <- spec$render$canvas$viewbox_width
      record$viewbox_height <- spec$render$canvas$viewbox_height
    } else {
      record$pixel_width <- journal$pixel_width
      record$pixel_height <- journal$pixel_height
      record$dpi <- journal$raster_dpi
      record$raster_source <- spec$generator_contract$outputs[[1L]]$path
    }
    record
  })
  identity$supporting_artifacts <- lapply(
    spec$generator_contract$supporting_artifacts,
    function(artifact) list(
      format = artifact$format,
      path = artifact$path,
      sha256 = workflow_sha256(staged$markdown)
    )
  )
  jsonlite::write_json(
    identity, staged$manifest, pretty = TRUE, auto_unbox = TRUE, null = "null"
  )
}

rendered_validation_arguments <- function(paths, staged, project_root) {
  c(
    paths$rendered_validator,
    "--project-root", project_root,
    "--spec", paths$spec,
    "--schema", paths$schema,
    "--generator", paths$generator,
    "--library", paths$library,
    "--rasterizer", paths$rasterizer,
    "--svg", staged$svg,
    "--png", staged$png,
    "--caption", staged$markdown,
    "--manifest", staged$manifest
  )
}

promote_artifacts <- function(spec, staged, project_root) {
  files <- c(
    stats::setNames(
      vapply(spec$generator_contract$outputs, function(item) item$path, character(1)),
      vapply(spec$generator_contract$outputs, function(item) item$format, character(1))
    ),
    markdown = spec$generator_contract$supporting_artifacts[[1L]]$path,
    manifest = spec$generator_contract$freshness_manifest
  )
  for (name in names(files)) {
    destination <- file.path(project_root, files[[name]])
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(staged[[name]], destination, overwrite = TRUE, copy.date = TRUE)) {
      stop(sprintf("Could not promote validated artifact: %s", destination), call. = FALSE)
    }
  }
}

main <- function() {
project_root <- find_project_root()
paths <- workflow_default_paths(project_root)
spec <- read_workflow_spec(paths$spec, paths$schema, validate_schema = TRUE)
validation <- validate_workflow_spec(spec, project_root, fail_on_render_blockers = TRUE)
assert_workflow_valid(validation)
if (length(spec$requires_author_decision) > 0L) {
  stop("Workflow rendering requires zero unresolved author decisions.", call. = FALSE)
}

node_runtime <- find_runtime("WORKFLOW_NODE", "node")
python_runtime <- find_runtime("WORKFLOW_PYTHON", "python")
sharp_module <- Sys.getenv("WORKFLOW_SHARP_MODULE", unset = "sharp")
journal <- workflow_index_by_id(spec$render$profiles)[["journal_landscape"]]

staging_parent <- file.path(paths$workflow_dir, ".staging")
dir.create(staging_parent, recursive = TRUE, showWarnings = FALSE)
staging_directory <- tempfile("workflow-render-", tmpdir = staging_parent)
dir.create(staging_directory, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(staging_directory, recursive = TRUE, force = TRUE), add = TRUE)

staged <- list(
  svg = file.path(staging_directory, basename(spec$generator_contract$outputs[[1L]]$path)),
  png = file.path(staging_directory, basename(spec$generator_contract$outputs[[2L]]$path)),
  markdown = file.path(staging_directory, basename(spec$generator_contract$supporting_artifacts[[1L]]$path)),
  manifest = file.path(staging_directory, basename(spec$generator_contract$freshness_manifest))
)

svg_lines <- svg_document(spec)
caption_lines <- caption_document(spec)
write_utf8_lines(svg_lines, staged$svg)
write_utf8_lines(caption_lines, staged$markdown)
message("Staged SVG SHA-256: ", workflow_sha256(staged$svg))
run_checked(
  node_runtime,
  c(
    paths$rasterizer, staged$svg, staged$png,
    as.character(journal$pixel_width), as.character(journal$pixel_height),
    as.character(journal$raster_dpi), sharp_module
  ),
  "Canonical SVG rasterization"
)
write_staged_manifest(spec, paths, staged)
run_checked(
  python_runtime,
  rendered_validation_arguments(paths, staged, project_root),
  "Staged workflow validation"
)

promote_artifacts(spec, staged, project_root)

canonical <- list(
  svg = file.path(project_root, spec$generator_contract$outputs[[1L]]$path),
  png = file.path(project_root, spec$generator_contract$outputs[[2L]]$path),
  markdown = file.path(project_root, spec$generator_contract$supporting_artifacts[[1L]]$path),
  manifest = file.path(project_root, spec$generator_contract$freshness_manifest)
)
run_checked(
  python_runtime,
  rendered_validation_arguments(paths, canonical, project_root),
  "Promoted workflow validation"
)
freshness <- validate_workflow_freshness(spec, paths)
if (!isTRUE(freshness$valid)) {
  stop(paste(c("Freshness validation failed:", freshness$errors), collapse = "\n"), call. = FALSE)
}

message("Canonical workflow assets promoted after all checks passed.")
for (name in names(canonical)) {
  message(sprintf("%s  %s", workflow_sha256(canonical[[name]]), canonical[[name]]))
}
}

main()
