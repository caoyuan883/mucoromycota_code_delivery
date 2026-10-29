#!/usr/bin/env Rscript
# BUSCO G/R/T comparison for 12 strains.
# Styling adapted from the supplied A4 v11 script: Times, horizontal bars,
# identical S/D/F/M colors, fixed bar height and one A4 page.
# Uses base R + grid. No ggplot2, Cairo or XQuartz needed for the PDF.
# In RStudio: open this script, click Source, then choose the input TSV.

## ---- [release adaptation] command-line options --------------------------------
## Added for the code release (see figures/plot_sources.tsv). Without options the
## original interactive behaviour is unchanged. Data handling and styling are unchanged.
.release_usage <- "Usage: Rscript plot_busco_genome_reference_study_12.R --input=TSV [--outdir=PARENT_DIR]\nTSV: busco_12strains_GTR_long_v1.tsv with columns sample, sample_order, subphylum, mode, lineage,\n     dataset_label, metric, percent, complete_percent, n_buscos\nDefault parent folder: folder of TSV; a new BUSCO_12_GRT_<timestamp> folder is created inside it.\nOutputs: Fig_BUSCO_12_GRT_A4.pdf/.png, check CSVs, figure_caption.txt, R_sessionInfo.txt"
.release_args <- commandArgs(trailingOnly = TRUE)
.release_opt <- function(name, default = NULL) {
  key <- paste0("--", name, "=")
  hit <- .release_args[startsWith(.release_args, key)]
  if (length(hit)) return(path.expand(sub(key, "", hit[length(hit)], fixed = TRUE)))
  default
}
if (any(.release_args %in% c("--help", "-h"))) {
  cat(.release_usage, "\n", sep = "")
  quit(save = "no", status = 0)
}
## ---- [end release adaptation] ---------------------------------------------------

# ---------------- Editable settings ----------------
INPUT_FILE <- .release_opt("input", "")                   # [release] --input; Empty: choose a file interactively.
OUTPUT_PARENT <- .release_opt("outdir", "")                # [release] --outdir; Empty: save beside the input table.
PAGE_W_MM <- 210
PAGE_H_MM <- 297
LEFT_MM <- 15
TOP_MM <- 13
BAR_SLOT_MM <- 3.4
BAR_FILL_FRACTION <- 0.82
STRAIN_GAP_MM <- 1.2
PANEL_GAP_MM <- 4
NAME_W_MM <- 52
MODE_W_MM <- 8
CODE_W_MM <- 6
PANEL_W_MM <- 103
C_GAP_MM <- 2
C_W_MM <- 9
FONT_FAMILY <- "Times"
TITLE_PT <- 11
SUBTITLE_PT <- 8
PANEL_TITLE_PT <- 9
SPECIES_PT <- 7.5
AXIS_PT <- 7
BAR_LABEL_PT <- 6
CAPTION_PT <- 8
LABEL_MIN <- 5
SHOW_COMPLETE_COLUMN <- TRUE
SHOW_CAPTION <- TRUE
SHOW_IN_PLOTS <- TRUE
SAVE_PNG <- TRUE
SAVE_SVG <- FALSE                  # Optional; base svg() needs Cairo support.
PNG_DPI <- 400
ROUNDING_TOLERANCE <- 0.21          # Percentage points, not a rescaling rule.
FIGURE_LABEL <- "Fig. X."
SUBPHYLA <- c("Mucoromycotina", "Mortierellomycotina", "Glomeromycotina")
DATASETS <- c("G", "R", "T")
METRICS <- c("S", "D", "F", "M")
COLORS <- c(S = "#2B8CBE", D = "#084081", F = "#FDAE61", M = "#D7191C")
CATEGORY_LABELS <- c(S = "Complete (Single-copy)", D = "Complete (Duplicated)",
                     F = "Fragmented", M = "Missing")

suppressPackageStartupMessages(library(grid))
if (INPUT_FILE == "") {
  if (!interactive()) stop("Set INPUT_FILE before running with Rscript.")
  message("Choose busco_12strains_GTR_long_v1.tsv")
  INPUT_FILE <- file.choose()
}
INPUT_FILE <- normalizePath(INPUT_FILE, mustWork = TRUE)
ext <- tolower(tools::file_ext(INPUT_FILE))
if (ext %in% c("xlsx", "xls")) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop('For Excel input, first run install.packages("readxl"); TSV requires no package.')
  raw <- as.data.frame(readxl::read_excel(INPUT_FILE, sheet = 1))
} else if (ext == "csv") {
  raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE,
                  fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
} else {
  raw <- read.delim(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"),
                    quote = "", comment.char = "")
}
required <- c("sample", "sample_order", "subphylum", "mode", "lineage",
              "dataset_label", "metric", "percent", "complete_percent", "n_buscos")
missing_cols <- setdiff(required, names(raw))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
d <- raw[, required]
for (nm in c("sample", "subphylum", "mode", "lineage", "dataset_label", "metric")) {
  d[[nm]] <- trimws(as.character(d[[nm]]))
  if (anyNA(d[[nm]]) || any(d[[nm]] == "")) stop("Missing values in ", nm)
}
d$metric <- toupper(d$metric)
d$dataset_label <- toupper(d$dataset_label)
for (nm in c("sample_order", "percent", "complete_percent", "n_buscos")) {
  d[[nm]] <- suppressWarnings(as.numeric(gsub("%", "", as.character(d[[nm]]), fixed = TRUE)))
  if (any(!is.finite(d[[nm]]))) stop("Invalid numeric values in ", nm)
}
if (!all(d$subphylum %in% SUBPHYLA)) stop("Unknown subphylum.")
if (!all(d$dataset_label %in% DATASETS)) stop("dataset_label must contain G, R and T.")
if (!all(d$metric %in% METRICS)) stop("metric must contain S, D, F and M.")
if (!all(d$mode %in% c("B2", "B3"))) stop("mode must be B2 or B3.")
if (any(d$percent < 0 | d$percent > 100 | d$complete_percent < 0 | d$complete_percent > 100))
  stop("Percentages must be in the range 0 to 100.")
if (!all(d$n_buscos == 1560)) stop("This figure expects n_buscos = 1560 for every row.")
if (!all(d$lineage %in% c("mucoromycota", "mucoromycota_odb12")))
  stop("This figure expects mucoromycota_odb12.")
if (any(d$sample_order <= 0 | d$sample_order != floor(d$sample_order)))
  stop("sample_order must be a positive integer.")
if (anyDuplicated(d[, c("sample", "dataset_label", "metric")]))
  stop("Duplicate sample/dataset/metric rows: resolve them; they will not be silently averaged.")
info <- unique(d[, c("sample", "sample_order", "subphylum", "mode")])
if (anyDuplicated(info$sample)) stop("Inconsistent sample metadata.")
if (nrow(info) != 12L) stop("Expected 12 strains; found ", nrow(info), ".")
if (anyDuplicated(info$sample_order)) stop("sample_order must be unique across strains.")
info <- info[order(match(info$subphylum, SUBPHYLA), info$sample_order), ]
rownames(info) <- NULL
counts <- table(factor(info$subphylum, levels = SUBPHYLA))
if (!identical(as.integer(counts), c(6L, 4L, 2L))) stop("Expected subphylum counts of 6, 4 and 2.")

bars_list <- list()
for (i in seq_len(nrow(info))) {
  for (ds in DATASETS) {
    z <- d[d$sample == info$sample[i] & d$dataset_label == ds, ]
    if (nrow(z) != 4L || !setequal(z$metric, METRICS))
      stop("Incomplete S/D/F/M data for ", info$sample[i], " / ", ds)
    if (length(unique(z$complete_percent)) != 1L)
      stop("Inconsistent complete_percent within a bar: ", info$sample[i], " / ", ds)
    z <- z[match(METRICS, z$metric), ]
    bars_list[[length(bars_list) + 1L]] <- data.frame(
      sample = info$sample[i], sample_order = info$sample_order[i],
      subphylum = info$subphylum[i], mode = info$mode[i], dataset_label = ds,
      S = z$percent[1], D = z$percent[2], F = z$percent[3], M = z$percent[4],
      C_reported = z$complete_percent[1], stringsAsFactors = FALSE)
  }
}
bars <- do.call(rbind, bars_list)
bars$total_SDFM <- rowSums(bars[, METRICS])
bars$total_minus_100 <- bars$total_SDFM - 100
bars$SD_minus_C <- bars$S + bars$D - bars$C_reported
bars$rounding_flag <- abs(bars$total_minus_100) > 1e-7 | abs(bars$SD_minus_C) > 1e-7
bars$check <- ifelse(abs(bars$total_minus_100) > ROUNDING_TOLERANCE |
                    abs(bars$SD_minus_C) > ROUNDING_TOLERANCE, "REVIEW", "PASS")
parent <- if (OUTPUT_PARENT == "") dirname(INPUT_FILE) else OUTPUT_PARENT
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
base_dir <- file.path(parent, paste0("BUSCO_12_GRT_", format(Sys.time(), "%Y%m%d_%H%M%S")))
OUT_DIR <- base_dir
j <- 1L
while (dir.exists(OUT_DIR)) {
  OUT_DIR <- paste0(base_dir, "_", j)
  j <- j + 1L
}
dir.create(OUT_DIR, recursive = TRUE)
write.csv(bars, file.path(OUT_DIR, "BUSCO_bar_data_check.csv"), row.names = FALSE)
if (any(bars$check == "REVIEW"))
  stop("Inconsistent percentages beyond rounding tolerance. See BUSCO_bar_data_check.csv.")
if (any(bars$rounding_flag))
  message("Small rounding differences found. Original category values and reported C are retained; no renormalization.")
X_MAX <- max(100.2, max(bars$total_SDFM))

wide <- info
for (ds in DATASETS) {
  z <- bars[bars$dataset_label == ds, ]
  wide[[ds]] <- z$C_reported[match(wide$sample, z$sample)]
}
wide$T_minus_R_pp <- round(wide$T - wide$R, 10)
wide$T_minus_G_pp <- round(wide$T - wide$G, 10)
wide$R_minus_G_pp <- round(wide$R - wide$G, 10)
write.csv(wide, file.path(OUT_DIR, "BUSCO_12_GRT_paired_summary.csv"), row.names = FALSE)
summary_names <- c("G", "R", "T", "T_minus_R_pp", "T_minus_G_pp", "R_minus_G_pp")
summary_table <- do.call(rbind, lapply(summary_names, function(nm) {
  x <- wide[[nm]]
  data.frame(metric = nm, n = length(x), median = median(x),
             minimum = min(x), maximum = max(x), stringsAsFactors = FALSE)
}))
write.csv(summary_table, file.path(OUT_DIR, "BUSCO_summary_statistics.csv"), row.names = FALSE)
write.csv(data.frame(
  comparison = c("T_gt_R", "T_eq_R", "T_lt_R", "abs_T_minus_G_le_1pp"),
  n_strains = c(sum(wide$T_minus_R_pp > 1e-7), sum(abs(wide$T_minus_R_pp) <= 1e-7),
                sum(wide$T_minus_R_pp < -1e-7), sum(abs(wide$T_minus_G_pp) <= 1 + 1e-7))),
  file.path(OUT_DIR, "BUSCO_paired_counts.csv"), row.names = FALSE)

# Absolute A4 layout; all three panels have the same width and bar height.
px <- LEFT_MM + NAME_W_MM + MODE_W_MM + CODE_W_MM
cx <- px + PANEL_W_MM + C_GAP_MM + C_W_MM / 2
mx <- LEFT_MM + NAME_W_MM + MODE_W_MM / 2
dx <- px - CODE_W_MM / 2
content_width <- NAME_W_MM + MODE_W_MM + CODE_W_MM + PANEL_W_MM + C_GAP_MM + C_W_MM
if (LEFT_MM + content_width > PAGE_W_MM - 10) stop("Configured figure width exceeds the page.")
plot_top <- TOP_MM + 26
layout_panels <- list()
y_cursor <- plot_top
for (sp in SUBPHYLA) {
  ix <- which(info$subphylum == sp)
  n <- length(ix)
  h <- n * 3 * BAR_SLOT_MM + (n - 1) * STRAIN_GAP_MM
  layout_panels[[sp]] <- list(indices = ix, heading = y_cursor,
                             top = y_cursor + 6, height = h)
  y_cursor <- y_cursor + 6 + h + 8 + PANEL_GAP_MM
}
axis_title_y <- y_cursor - PANEL_GAP_MM + 1
caption_top <- axis_title_y + 8
if (caption_top + ifelse(SHOW_CAPTION, 25, 0) > PAGE_H_MM - 10)
  stop("Figure exceeds A4 height. Reduce BAR_SLOT_MM or CAPTION_PT.")

caption <- paste0(
  FIGURE_LABEL, " BUSCO completeness of 12 paired genome assemblies and proteomes. ",
  "G, genome; R, public reference proteome; T, this study's delivered representative proteome. ",
  "All datasets use mucoromycota_odb12 (n = 1,560). Bars retain the reported S, D, F and M percentages. ",
  "B2/B3 identifies the annotation mode used in this study. ",
  if (SHOW_COMPLETE_COLUMN) "The right-hand column gives reported complete BUSCO values (C). " else "",
  "Category totals or S + D may differ slightly from 100% or reported C because values are rounded. ",
  sprintf("Segments of at least %g%% are labelled.", LABEL_MIN))
writeLines(caption, file.path(OUT_DIR, "figure_caption.txt"), useBytes = TRUE)

plot_busco <- function() {
  grid.newpage()
  dw <- convertWidth(unit(1, "npc"), "mm", valueOnly = TRUE)
  dh <- convertHeight(unit(1, "npc"), "mm", valueOnly = TRUE)
  sc <- min(dw / PAGE_W_MM, dh / PAGE_H_MM)
  pushViewport(viewport(width = unit(PAGE_W_MM * sc, "mm"),
                        height = unit(PAGE_H_MM * sc, "mm")))
  on.exit(popViewport(), add = TRUE)
  xx <- function(x) unit(x / PAGE_W_MM, "npc")
  yy <- function(y) unit(1 - y / PAGE_H_MM, "npc")
  txt <- function(label, x, y, size = AXIS_PT, face = "plain", col = "black", just = "centre") {
    grid.text(label, x = xx(x), y = yy(y), just = just,
              gp = gpar(fontfamily = FONT_FAMILY, fontsize = size * sc,
                        fontface = face, col = col))
  }
  rect <- function(x, y, w, h, fill, border = NA, lwd = 0.6) {
    grid.rect(x = xx(x + w / 2), y = yy(y + h / 2),
              width = unit(w / PAGE_W_MM, "npc"), height = unit(h / PAGE_H_MM, "npc"),
              gp = gpar(fill = fill, col = border, lwd = lwd * sc))
  }
  line <- function(x0, y0, x1, y1, col = "black", lwd = 0.5, lty = "solid") {
    grid.lines(x = xx(c(x0, x1)), y = yy(c(y0, y1)),
               gp = gpar(col = col, lwd = lwd * sc, lty = lty))
  }
  tx_width <- function(label, size, face = "plain") {
    g <- textGrob(label, gp = gpar(fontfamily = FONT_FAMILY, fontsize = size * sc, fontface = face))
    convertWidth(grobWidth(g), "mm", valueOnly = TRUE) / sc
  }
  txt("BUSCO completeness of 12 genome-reference-study triplets",
      PAGE_W_MM / 2, TOP_MM, TITLE_PT, "bold")
  txt("mucoromycota_odb12 (n = 1,560)", PAGE_W_MM / 2, TOP_MM + 5.5, SUBTITLE_PT)
  widths <- vapply(CATEGORY_LABELS, function(z) tx_width(z, 7), numeric(1))
  legend_w <- sum(widths + 8) + 3 * 3
  lx <- LEFT_MM + (content_width - legend_w) / 2
  if (legend_w > content_width) stop("Legend is too wide; reduce legend font size.")
  for (k in seq_along(METRICS)) {
    rect(lx, TOP_MM + 10.5, 5, 2.5, COLORS[[METRICS[k]]])
    txt(CATEGORY_LABELS[k], lx + 6.5, TOP_MM + 11.75, 7, just = "left")
    lx <- lx + 8 + widths[k] + 3
  }
  txt("G = genome   R = public reference   T = this study   |   B2/B3 = study mode",
      PAGE_W_MM / 2, TOP_MM + 19, 7.5)

  for (p in seq_along(SUBPHYLA)) {
    sp <- SUBPHYLA[p]
    lp <- layout_panels[[sp]]
    txt(sprintf("%s  %s (n = %d)", letters[p], sp, length(lp$indices)),
        LEFT_MM, lp$heading, PANEL_TITLE_PT, "bold", just = "left")
    txt("Mode", mx, lp$heading, 7, "bold")
    if (SHOW_COMPLETE_COLUMN) txt("C (%)", cx, lp$heading, 7, "bold")
    for (tick in seq(0, 100, 20)) {
      xt <- px + tick / X_MAX * PANEL_W_MM
      line(xt, lp$top, xt, lp$top + lp$height, "grey85", 0.5, "dashed")
    }
    for (j in seq_along(lp$indices)) {
      i <- lp$indices[j]
      sname <- info$sample[i]
      y0 <- lp$top + (j - 1) * (3 * BAR_SLOT_MM + STRAIN_GAP_MM)
      display_name <- gsub("_", " ", sname, fixed = TRUE)
      sz <- SPECIES_PT
      while (tx_width(display_name, sz, "italic") > NAME_W_MM - 1 && sz > 6)
        sz <- sz - 0.1
      if (tx_width(display_name, sz, "italic") > NAME_W_MM - 1)
        stop("Species name does not fit. Increase NAME_W_MM: ", display_name)
      txt(display_name, LEFT_MM, y0 + 1.5 * BAR_SLOT_MM, sz, "italic", just = "left")
      txt(info$mode[i], mx, y0 + 1.5 * BAR_SLOT_MM, 7)
      for (k in seq_along(DATASETS)) {
        ds <- DATASETS[k]
        z <- bars[bars$sample == sname & bars$dataset_label == ds, ]
        v <- as.numeric(unlist(z[1, METRICS], use.names = FALSE))
        ends <- cumsum(v)
        starts <- c(0, head(ends, -1))
        yc <- y0 + (k - 0.5) * BAR_SLOT_MM
        bh <- BAR_SLOT_MM * BAR_FILL_FRACTION
        txt(ds, dx, yc, 7)
        for (m in seq_along(METRICS)) {
          if (v[m] > 0) rect(px + starts[m] / X_MAX * PANEL_W_MM, yc - bh / 2,
                             v[m] / X_MAX * PANEL_W_MM, bh, COLORS[[METRICS[m]]])
          if (v[m] >= LABEL_MIN)
            txt(sprintf("%.1f", v[m]), px + (starts[m] + v[m] / 2) / X_MAX * PANEL_W_MM,
                yc, BAR_LABEL_PT, "bold", if (METRICS[m] == "D") "white" else "black")
        }
        if (SHOW_COMPLETE_COLUMN) txt(sprintf("%.1f", z$C_reported), cx, yc, 7)
      }
    }
    rect(px, lp$top, PANEL_W_MM, lp$height, NA, "black", 0.7)
    for (tick in seq(0, 100, 20)) {
      xt <- px + tick / X_MAX * PANEL_W_MM
      line(xt, lp$top + lp$height, xt, lp$top + lp$height + 1, lwd = 0.5)
      txt(as.character(tick), xt, lp$top + lp$height + 3.6, AXIS_PT)
    }
  }
  txt("BUSCO percentage (%)", px + PANEL_W_MM / 2, axis_title_y, 8, "bold")
  if (SHOW_CAPTION) {
    words <- strsplit(caption, " +")[[1]]
    cap_lines <- character()
    current <- ""
    for (w in words) {
      candidate <- if (current == "") w else paste(current, w)
      if (tx_width(candidate, CAPTION_PT) > content_width && current != "") {
        cap_lines <- c(cap_lines, current)
        current <- w
      } else current <- candidate
    }
    cap_lines <- c(cap_lines, current)
    lh <- CAPTION_PT * 25.4 / 72 * 1.22
    if (caption_top + length(cap_lines) * lh > PAGE_H_MM - 10)
      stop("Caption does not fit; reduce CAPTION_PT or BAR_SLOT_MM.")
    for (i in seq_along(cap_lines))
      txt(cap_lines[i], LEFT_MM, caption_top + (i - 1) * lh, CAPTION_PT, just = "left")
  }
}

save_graphic <- function(path, kind) {
  old <- dev.cur()
  success <- FALSE
  if (kind == "pdf") {
    grDevices::pdf(path, width = PAGE_W_MM / 25.4, height = PAGE_H_MM / 25.4,
                   family = FONT_FAMILY, useDingbats = FALSE)
  } else if (kind == "png") {
    grDevices::png(path, width = PAGE_W_MM, height = PAGE_H_MM, units = "mm", res = PNG_DPI)
  } else {
    grDevices::svg(path, width = PAGE_W_MM / 25.4, height = PAGE_H_MM / 25.4,
                   family = FONT_FAMILY)
  }
  this_device <- dev.cur()
  on.exit({
    if (this_device %in% dev.list()) dev.off(this_device)
    if (!success && file.exists(path)) unlink(path)
  }, add = TRUE)
  plot_busco()
  dev.off(this_device)
  success <- file.exists(path) && file.info(path)$size > 0
  if (!success) stop("Output was not created: ", path)
  message("Saved: ", path)
}
PDF_FILE <- file.path(OUT_DIR, "Fig_BUSCO_12_GRT_A4.pdf")
save_graphic(PDF_FILE, "pdf")
if (SAVE_PNG) {
  tryCatch(save_graphic(file.path(OUT_DIR, "Fig_BUSCO_12_GRT_A4.png"), "png"),
           error = function(e) warning("PNG export failed; PDF retained: ", conditionMessage(e)))
}
if (SAVE_SVG) {
  if (capabilities("cairo")) {
    tryCatch(save_graphic(file.path(OUT_DIR, "Fig_BUSCO_12_GRT_A4.svg"), "svg"),
             error = function(e) warning("SVG export failed; PDF retained: ", conditionMessage(e)))
  } else message("SVG skipped: Cairo is unavailable. The base-R PDF does not require Cairo.")
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "R_sessionInfo.txt"))
writeLines(c(paste("Input:", INPUT_FILE), paste("Input MD5:", unname(tools::md5sum(INPUT_FILE))),
             "C statistics use complete_percent, not S + D recalculated from rounded categories.",
             "Bar widths retain original percentages; categories are not normalized to sum to 100.",
             "Statistics are across 12 strain-level values, not pooled across proteins.",
             "Within each subphylum, sample_order controls the display order."),
           file.path(OUT_DIR, "analysis_notes.txt"))
message("\nSummary based on reported C:")
print(summary_table, row.names = FALSE)
message("Output directory: ", normalizePath(OUT_DIR))
if (SHOW_IN_PLOTS && interactive()) {
  tryCatch(plot_busco(), error = function(e) message("Preview unavailable; open the exported PDF: ", conditionMessage(e)))
}
