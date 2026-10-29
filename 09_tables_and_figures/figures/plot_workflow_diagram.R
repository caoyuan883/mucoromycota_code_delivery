#!/usr/bin/env Rscript
# Mucoromycota structural annotation resource - workflow scheme (Methods), v4
# v4: Helvetica Bold throughout; robust PDF export on macOS
#     (cairo_pdf -> quartz pdf -> base pdf fallback, file existence checked).
# Canvas: 140 x 172 mm; Helvetica Bold 7 pt (Scientific Data)

## ---- [release adaptation] command-line options --------------------------------
## Added for the code release (see figures/plot_sources.tsv). Without options the
## original interactive behaviour is unchanged. Data handling and styling are unchanged.
.release_usage <- "Usage: Rscript plot_workflow_diagram.R [--outdir=DIR]\nNo input file: boxes and text are defined inside the script.\nDefault output folder: RStudio folder chooser if available, else current working directory.\nOutputs: Fig_workflow_v4.pdf, Fig_workflow_v4.png, Fig_workflow_v4.svg (svg only if svglite is installed)"
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

packages <- c("ggplot2", "tibble")
miss <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) stop("Missing R packages (install them first; this script does not install packages): ", paste(miss, collapse = ", "), call. = FALSE)  # [release] was: install.packages()
library(ggplot2); library(tibble)

out_dir <- if (!is.null(.release_opt("outdir"))) {  # [release] --outdir first
  dir.create(.release_opt("outdir"), showWarnings = FALSE, recursive = TRUE); .release_opt("outdir")
} else if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::selectDirectory(caption = "Select output folder for the workflow figure")
} else getwd()
if (is.null(out_dir) || !dir.exists(out_dir)) out_dir <- getwd()

W <- 140; H <- 172
fs <- 7.0; font <- "Helvetica"; face <- "bold"
line_col <- "#333333"; dec_fill <- "#FDD9B5"; frame_col <- "#9A9A9A"
lw <- 0.3
arr <- arrow(length = unit(1.6, "mm"), type = "closed", angle = 22)

# ---------------- nodes ----------------
nodes <- tribble(
  ~id,      ~cx, ~cy,  ~w,  ~h, ~kind, ~label,
  "ncbi",    40, 164,  38,   8, "box", "NCBI genomes",
  "jgi",    100, 164,  42,   8, "box", "JGI MycoCosm genomes",
  "cand",    70, 147,  38,   8, "box", "candidate genomes",
  "audit",   70, 133,  46,  10, "box", "taxonomy, redundancy\nand BUSCO audit",
  "cohort",  70, 119,  38,   8, "box", "study cohort",
  "rep",     23, 101,  40,  10, "box", "RepeatModeler and\nRepeatMasker",
  "sra",     70, 101,  38,  10, "box", "RNA-seq from SRA",
  "prot",   117, 101,  40,  10, "box", "OrthoDB v12\nfungal proteins",
  "mask",    23,  86,  40,   8, "box", "soft-masked genome",
  "rna",     70,  86,  38,  10, "box", "align and filter\nRNA-seq libraries",
  "dec",     70,  70,  42,  14, "dec", "RNA-seq\nretained?",
  "b3",      42,  46,  32,   8, "box", "run BRAKER3",
  "b2",      98,  46,  32,   8, "box", "run BRAKER2",
  "ann",     70,  25,  52,  10, "box", "gene annotation\n(longest isoform per gene)",
  "qc",      22,   7,  38,  10, "box", "quality control\n(BUSCO, OMArk)",
  "fun",     70,   7,  52,  10, "box", "functional annotation\n(InterProScan, eggNOG-mapper)",
  "phy",    118,   7,  38,  10, "box", "phylogenomic\nvalidation"
)

frame <- c(x0 = 5, x1 = 135, y0 = 35, y1 = 58)

# ---------------- connectors ----------------
lines <- tribble(
  ~x, ~y, ~xend, ~yend,
  40, 160,  40, 156,
  100, 160, 100, 156,
  40, 156, 100, 156,
  23, 111, 117, 111,
  49,  70,  42,  70,
  91,  70,  98,  70,
  70,  20,  70,  16,
  22,  16, 118,  16
)
arrows <- tribble(
  ~x, ~y, ~xend, ~yend, ~style,
  70, 156,  70, 151, "solid",
  70, 143,  70, 138, "solid",
  70, 128,  70, 123, "solid",
  70, 115,  70, 106, "solid",
  23, 111,  23, 106, "solid",
  117, 111, 117, 106, "solid",
  23,  96,  23,  90, "solid",
  70,  96,  70,  91, "solid",
  70,  81,  70,  77, "solid",
  42,  70,  42,  50, "solid",
  98,  70,  98,  50, "solid",
  23,  82,  23, frame[["y1"]], "solid",
  117, 96, 117, frame[["y1"]], "solid",
  58,  46,  82,  46, "dashed",
  70, frame[["y0"]], 70, 30, "solid",
  22,  16,  22,  12, "solid",
  70,  16,  70,  12, "solid",
  118, 16, 118,  12, "solid"
)

labels <- tribble(
  ~x,   ~y,  ~label,      ~fontface,
  45.5, 72.4, "yes",       "bold",
  94.5, 72.4, "no",        "bold",
  70,  48.4, "low BUSCO", "bold",
  70,  38.2, "gene prediction (BRAKER, fungal mode, soft-masked genome)", "bold.italic"
)

# ---------------- shapes ----------------
round_rect <- function(x0, y0, w, h, r = 1.2, n = 8) {
  r <- min(r, w / 2, h / 2)
  a <- function(cx, cy, from, to) {
    t <- seq(from, to, length.out = n)
    data.frame(x = cx + r * cos(t), y = cy + r * sin(t))
  }
  rbind(a(x0 + w - r, y0 + h - r, 0, pi / 2), a(x0 + r, y0 + h - r, pi / 2, pi),
        a(x0 + r, y0 + r, pi, 3 * pi / 2), a(x0 + w - r, y0 + r, 3 * pi / 2, 2 * pi))
}
diamond <- function(cx, cy, w, h)
  data.frame(x = c(cx, cx + w / 2, cx, cx - w / 2), y = c(cy + h / 2, cy, cy - h / 2, cy))

polys <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
  nd <- nodes[i, ]
  d <- if (nd$kind == "dec") diamond(nd$cx, nd$cy, nd$w, nd$h)
       else round_rect(nd$cx - nd$w / 2, nd$cy - nd$h / 2, nd$w, nd$h)
  d$id <- nd$id; d$fill <- ifelse(nd$kind == "dec", dec_fill, "white"); d
}))
frame_poly <- round_rect(frame[["x0"]], frame[["y0"]], frame[["x1"]] - frame[["x0"]],
                         frame[["y1"]] - frame[["y0"]], r = 2)

# ---------------- plot ----------------
p <- ggplot() +
  geom_polygon(data = frame_poly, aes(x, y), fill = NA, colour = frame_col,
               linewidth = lw, linetype = "33") +
  geom_segment(data = lines, aes(x = x, y = y, xend = xend, yend = yend),
               colour = line_col, linewidth = lw, lineend = "square") +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, linetype = style),
               colour = line_col, linewidth = lw, arrow = arr, linejoin = "mitre") +
  geom_polygon(data = polys, aes(x, y, group = id, fill = fill),
               colour = line_col, linewidth = lw) +
  geom_text(data = nodes, aes(cx, cy, label = label), size = fs / .pt,
            lineheight = 0.95, family = font, fontface = face) +
  geom_text(data = labels, aes(x, y, label = label, fontface = fontface),
            size = (fs - 0.5) / .pt, family = font, colour = "#222222") +
  scale_fill_identity() +
  scale_linetype_manual(values = c(solid = "solid", dashed = "22"), guide = "none") +
  coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") +
  theme_void(base_family = font) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(1, 1, 1, 1))

# ---------------- export ----------------
stem <- file.path(out_dir, "Fig_workflow_v4")
pdf_file <- paste0(stem, ".pdf"); png_file <- paste0(stem, ".png"); svg_file <- paste0(stem, ".svg")
if (file.exists(pdf_file)) file.remove(pdf_file)

save_pdf <- function() {
  try_dev <- function(dev, name) {
    ok <- tryCatch({
      suppressWarnings(ggsave(pdf_file, p, width = W, height = H, units = "mm", device = dev))
      file.exists(pdf_file) && file.info(pdf_file)$size > 0
    }, error = function(e) FALSE)
    if (ok) message("PDF written with ", name)
    ok
  }
  if (isTRUE(capabilities("cairo")) && try_dev(cairo_pdf, "cairo_pdf")) return(invisible(TRUE))
  if (Sys.info()[["sysname"]] == "Darwin" &&
      try_dev(function(filename, width, height, ...)
        grDevices::quartz(type = "pdf", file = filename, width = width, height = height),
        "quartz")) return(invisible(TRUE))
  if (try_dev(function(filename, width, height, ...)
        grDevices::pdf(filename, width = width, height = height, family = "Helvetica"),
        "base pdf")) return(invisible(TRUE))
  stop("PDF export failed with all devices")
}
save_pdf()

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(png_file, p, width = W, height = H, units = "mm", dpi = 600,
         device = ragg::agg_png, bg = "white")
} else {
  ggsave(png_file, p, width = W, height = H, units = "mm", dpi = 600, bg = "white")
}
if (requireNamespace("svglite", quietly = TRUE))
  ggsave(svg_file, p, width = W, height = H, units = "mm")

print(p)
message("PDF: ", pdf_file, "  exists = ", file.exists(pdf_file))
message("PNG: ", png_file)
