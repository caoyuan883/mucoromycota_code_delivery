# Fig. 1a: proportional database overlap and one annotation-availability bar
# Requires only ggplot2. No input file is required.
# Outputs are saved in the current working directory shown by getwd().

## ---- [release adaptation] command-line options --------------------------------
## Added for the code release (see figures/plot_sources.tsv). Without options the
## original interactive behaviour is unchanged. Data handling and styling are unchanged.
.release_usage <- "Usage: Rscript plot_database_overlap_annotation_availability.R [--outdir=DIR]\nNo input file: genome counts are constants inside the script (see figures/README.md).\nDefault output folder: current working directory.\nOutputs: Fig1a_FINAL_proportional_circles_single_bar.pdf/.tiff, *_preview.png, Fig1a_FINAL_counts.tsv"
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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is not installed.")
}
library(ggplot2)

output_dir <- .release_opt("outdir", getwd())  # [release] was: getwd()
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat("Outputs will be saved in:\n", output_dir, "\n\n", sep = "")

# Final verified counts.
n_ncbi_total <- 180
n_jgi_total <- 120
n_overlap <- 98
n_ncbi_only <- 82
n_jgi_only <- 22
n_selected <- 202

# Public protein availability:
# 120 JGI-associated genomes all have public proteins;
# 27 additional annotated genomes are NCBI-only;
# 55 NCBI-only genomes have no public protein set.
n_jgi_annotated <- 120
n_ncbi_only_annotated <- 27
n_unannotated <- 55
n_annotated_total <- n_jgi_annotated + n_ncbi_only_annotated

stopifnot(n_ncbi_only + n_overlap + n_jgi_only == n_selected)
stopifnot(n_jgi_annotated + n_ncbi_only_annotated + n_unannotated == n_selected)
stopifnot(n_annotated_total == 147)

circle_intersection_area <- function(r1, r2, d) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  a2 <- acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  term <- 0.5 * sqrt(
    (-d + r1 + r2) *
      (d + r1 - r2) *
      (d - r1 + r2) *
      (d + r1 + r2)
  )
  r1^2 * a1 + r2^2 * a2 - term
}

circle_polygon <- function(cx, cy, radius, group, database, n = 800) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = cx + radius * cos(theta),
    y = cy + radius * sin(theta),
    group = group,
    database = database,
    stringsAsFactors = FALSE
  )
}

# Circle area is proportional to database coverage.
area_scale <- 0.012
r_ncbi <- sqrt(n_ncbi_total * area_scale / pi)
r_jgi <- sqrt(n_jgi_total * area_scale / pi)
target_overlap_area <- n_overlap * area_scale

distance <- uniroot(
  function(d) circle_intersection_area(r_ncbi, r_jgi, d) - target_overlap_area,
  interval = c(abs(r_ncbi - r_jgi) + 1e-6, r_ncbi + r_jgi - 1e-6)
)$root

left_center <- -2.25 - distance / 2
right_center <- -2.25 + distance / 2
circle_y <- 0.15

circle_data <- rbind(
  circle_polygon(left_center, circle_y, r_ncbi, "NCBI", "NCBI"),
  circle_polygon(right_center, circle_y, r_jgi, "JGI", "JGI")
)

left_number_x <- left_center - 0.42 * r_ncbi
overlap_number_x <- (left_center + right_center) / 2
right_number_x <- right_center + 0.42 * r_jgi

# One bar representing all 202 selected genomes.
bar_left <- 0.55
bar_width <- 4.25
bar_ymin <- -0.18
bar_ymax <- 0.28

bar_segments <- data.frame(
  category = c(
    "JGI-associated public proteins",
    "NCBI-only public proteins",
    "No public protein set"
  ),
  value = c(n_jgi_annotated, n_ncbi_only_annotated, n_unannotated),
  colour = c("#E69F00", "#0072B2", "#D9D9D9"),
  stringsAsFactors = FALSE
)
bar_segments$xmin_count <- c(0, cumsum(bar_segments$value)[-nrow(bar_segments)])
bar_segments$xmax_count <- cumsum(bar_segments$value)
bar_segments$xmin <- bar_left + bar_width * bar_segments$xmin_count / n_selected
bar_segments$xmax <- bar_left + bar_width * bar_segments$xmax_count / n_selected
bar_segments$xmid <- (bar_segments$xmin + bar_segments$xmax) / 2
bar_segments$label <- c(
  "120 JGI-associated",
  "27 NCBI-only",
  "55 without"
)
bar_segments$text_colour <- c("#202020", "white", "#303030")

annotated_end <- bar_left + bar_width * n_annotated_total / n_selected
bar_right <- bar_left + bar_width

p <- ggplot() +
  geom_polygon(
    data = circle_data,
    aes(x = x, y = y, group = group, fill = database, colour = database),
    alpha = 0.20,
    linewidth = 1.0
  ) +
  annotate(
    "text", x = left_center - 0.18, y = circle_y + r_ncbi + 0.30,
    label = "NCBI", colour = "#0072B2",
    fontface = "bold", size = 5.0
  ) +
  annotate(
    "text", x = right_center + 0.16, y = circle_y + r_jgi + 0.30,
    label = "JGI", colour = "#D89000",
    fontface = "bold", size = 5.0
  ) +
  annotate("text", x = left_number_x, y = 0.20,
           label = "82", fontface = "bold", size = 7.0, colour = "#202020") +
  annotate("text", x = overlap_number_x, y = 0.20,
           label = "98", fontface = "bold", size = 7.0, colour = "#202020") +
  annotate("text", x = right_number_x, y = 0.20,
           label = "22", fontface = "bold", size = 7.0, colour = "#202020") +
  annotate("text", x = left_number_x, y = -0.12,
           label = "NCBI only", size = 3.0, colour = "#303030") +
  annotate("text", x = overlap_number_x, y = -0.12,
           label = "Both", size = 3.0, colour = "#303030") +
  annotate("text", x = right_number_x, y = -0.12,
           label = "JGI only", size = 3.0, colour = "#303030") +
  annotate(
    "text", x = -2.25, y = 1.67,
    label = "Database coverage",
    size = 4.5, fontface = "bold", colour = "#202020"
  ) +
  annotate(
    "text", x = -2.25, y = -1.20,
    label = "Selected assembly source: NCBI 180; JGI 22",
    size = 3.25, fontface = "bold", colour = "#404040"
  ) +
  annotate(
    "text", x = 2.68, y = 1.67,
    label = "Public protein availability",
    size = 4.5, fontface = "bold", colour = "#202020"
  ) +
  geom_rect(
    data = bar_segments,
    aes(xmin = xmin, xmax = xmax, ymin = bar_ymin, ymax = bar_ymax),
    inherit.aes = FALSE,
    fill = bar_segments$colour,
    colour = NA
  ) +
  annotate(
    "rect", xmin = bar_left, xmax = bar_right,
    ymin = bar_ymin, ymax = bar_ymax,
    fill = NA, colour = "#6F6F6F", linewidth = 0.45
  ) +
  geom_text(
    data = bar_segments,
    aes(x = xmid, y = (bar_ymin + bar_ymax) / 2, label = label),
    inherit.aes = FALSE,
    colour = bar_segments$text_colour,
    size = 2.75,
    fontface = "bold"
  ) +
  annotate(
    "segment", x = bar_left, xend = annotated_end,
    y = 0.64, yend = 0.64,
    colour = "#333333", linewidth = 0.42
  ) +
  annotate(
    "segment", x = bar_left, xend = bar_left,
    y = 0.56, yend = 0.72,
    colour = "#333333", linewidth = 0.42
  ) +
  annotate(
    "segment", x = annotated_end, xend = annotated_end,
    y = 0.56, yend = 0.72,
    colour = "#333333", linewidth = 0.42
  ) +
  annotate(
    "text", x = (bar_left + annotated_end) / 2, y = 0.86,
    label = "147 genomes with public protein sets",
    size = 3.45, fontface = "bold", colour = "#303030"
  ) +
  annotate(
    "text", x = 2.68, y = -0.60,
    label = "The 147 annotated genomes include all 120 JGI-associated genomes\nand 27 additional NCBI-only genomes.",
    size = 3.05, lineheight = 1.02, colour = "#404040"
  ) +
  annotate(
    "text", x = 2.68, y = -1.04,
    label = "The remaining 55 genomes are NCBI-only and lacked public protein annotations.",
    size = 3.05, colour = "#555555"
  ) +
  scale_fill_manual(values = c("NCBI" = "#56B4E9", "JGI" = "#E69F00")) +
  scale_colour_manual(values = c("NCBI" = "#0072B2", "JGI" = "#D89000")) +
  coord_fixed(xlim = c(-4.15, 5.25), ylim = c(-1.55, 1.95), clip = "off") +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(8, 10, 8, 10)
  )

pdf_file <- file.path(output_dir, "Fig1a_FINAL_proportional_circles_single_bar.pdf")
tiff_file <- file.path(output_dir, "Fig1a_FINAL_proportional_circles_single_bar.tiff")
png_file <- file.path(output_dir, "Fig1a_FINAL_proportional_circles_single_bar_preview.png")
tsv_file <- file.path(output_dir, "Fig1a_FINAL_counts.tsv")

ggsave(pdf_file, p, width = 7.2, height = 3.8, units = "in")
ggsave(
  tiff_file, p,
  width = 7.2, height = 3.8, units = "in",
  dpi = 600, compression = "lzw"
)
ggsave(png_file, p, width = 7.2, height = 3.8, units = "in", dpi = 300)

summary_table <- data.frame(
  category = c("NCBI only", "NCBI and JGI", "JGI only", "Overall"),
  total = c(82, 98, 22, 202),
  public_proteins = c(27, 98, 22, 147),
  no_public_proteins = c(55, 0, 0, 55),
  stringsAsFactors = FALSE
)
write.table(summary_table, tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat(
  "\nOutputs saved in:\n", output_dir, "\n",
  "Fig1a_FINAL_proportional_circles_single_bar.pdf\n",
  "Fig1a_FINAL_proportional_circles_single_bar.tiff\n",
  "Fig1a_FINAL_proportional_circles_single_bar_preview.png\n",
  "Fig1a_FINAL_counts.tsv\n",
  sep = ""
)

print(summary_table)
print(p)
