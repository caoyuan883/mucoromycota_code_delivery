#!/usr/bin/env Rscript
# Umbelopsis ramanniana BUSCO completeness figure (Methods, gene structure annotation) - v2
# Input : Umbelopsis_ramanniana_BUSCO_completeness-2.tsv (selected interactively)
# Output: <input>_plot_v2.pdf / .png / .svg + <input>_data_check.csv (same folder as input)
# Rows  : Genome, Public annotation, BRAKER3 (longest isoform), BRAKER2 (longest isoform, final)
#         BRAKER2 raw braker.aa row (compleasm 0.2.8) is excluded; only BUSCO v6.0.0 rows plotted.
# Colour: same palette as the main longest-isoform BUSCO figure (blue / navy / orange / red).
#         Checked by CVD simulation (Machado 2009, severity 1.0): min pairwise dE76 = 27.8
#         (deuteranopia); CIELAB L* = 55.1 / 27.6 / 77.3 / 45.9, so categories stay separable in greyscale.

packages <- c("ggplot2", "dplyr", "tidyr", "readr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) install.packages(missing_packages, repos = "https://cloud.r-project.org")
library(ggplot2); library(dplyr); library(tidyr); library(readr)

message("Select Umbelopsis_ramanniana_BUSCO_completeness*.tsv")
input_file <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::selectFile(caption = "Select Umbelopsis BUSCO TSV", filter = "TSV (*.tsv)")
} else file.choose()
stopifnot(!is.null(input_file), file.exists(input_file))
out_dir <- dirname(input_file)
input_stem <- tools::file_path_sans_ext(basename(input_file))

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
dat$dataset <- trimws(dat$dataset); dat$tool <- trimws(dat$tool)

dat <- dat %>%
  filter(grepl("^BUSCO", tool), !grepl("raw", dataset, ignore.case = TRUE)) %>%
  mutate(label = case_when(
    grepl("^Genome",  dataset) ~ "Genome (G)",
    grepl("^Public",  dataset) ~ "Public annotation (P)",
    grepl("^BRAKER3", dataset) ~ "BRAKER3 protein (B3)",
    grepl("^BRAKER2", dataset) ~ "BRAKER2 protein (B2)*"
  ))
stopifnot(nrow(dat) == 4, !any(is.na(dat$label)))

check <- dat %>%
  transmute(label, C_pct, S_plus_D = S_pct + D_pct, total = S_pct + D_pct + F_pct + M_pct)
write_csv(check, file.path(out_dir, paste0(input_stem, "_data_check.csv")))
print(check)
if (any(check$total < 99.5 | check$total > 100.5)) warning("S+D+F+M not close to 100%; see _data_check.csv")

y_order <- c("Genome (G)", "Public annotation (P)", "BRAKER3 protein (B3)", "BRAKER2 protein (B2)*")
legend_order <- c("Complete (Single-copy)", "Complete (Duplicated)", "Fragmented", "Missing")

colors <- c(
  "Complete (Single-copy)" = "#2B8CBE",
  "Complete (Duplicated)"  = "#084081",
  "Fragmented"             = "#FDAE61",
  "Missing"                = "#D7191C"
)
text_colors <- c(
  "Complete (Single-copy)" = "white",
  "Complete (Duplicated)"  = "white",
  "Fragmented"             = "black",
  "Missing"                = "white"
)

long <- dat %>%
  pivot_longer(c(S_pct, D_pct, F_pct, M_pct), names_to = "metric", values_to = "percent") %>%
  mutate(
    category = recode(metric, S_pct = legend_order[1], D_pct = legend_order[2],
                      F_pct = legend_order[3], M_pct = legend_order[4]),
    category = factor(category, levels = legend_order),
    label = factor(label, levels = rev(y_order)),
    txt = ifelse(percent >= 5, sprintf("%.1f", percent), ""),
    txt_col = unname(text_colors[as.character(category)])
  )

cvals <- dat %>%
  mutate(label = factor(label, levels = rev(y_order)), ctxt = sprintf("C = %.1f%%", C_pct))

n_busco <- unique(dat$n_buscos)
lineage <- unique(dat$lineage)

p <- ggplot(long, aes(x = percent, y = label, fill = category)) +
  geom_col(width = 0.72, colour = NA, position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = txt, colour = txt_col),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3.0, fontface = "bold", family = "serif", show.legend = FALSE) +
  geom_text(data = cvals, aes(x = 101, y = label, label = ctxt), inherit.aes = FALSE,
            hjust = 0, size = 3.0, fontface = "bold", family = "serif") +
  scale_colour_identity() +
  scale_fill_manual(values = colors, breaks = legend_order, drop = FALSE) +
  scale_x_continuous(limits = c(0, 116), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = expression(bolditalic("Umbelopsis ramanniana")),
    x = "BUSCO percentage (%)", y = NULL, fill = "BUSCO category",
    caption = paste0("BUSCO v6.0.0, ", lineage, " (n = ", n_busco, "). ",
                     "B3 and B2 protein sets filtered to the longest isoform per gene. ",
                     "*Final annotation.")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 10, base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12.5),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 9, colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.35, linetype = "dashed"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(1.1, "cm"),
    plot.caption = element_text(hjust = 0, size = 7.5),
    plot.margin = margin(8, 12, 8, 8)
  )

fig_width <- 7.2
fig_height <- 3.4

pdf_file <- file.path(out_dir, paste0(input_stem, "_plot_v2.pdf"))
png_file <- file.path(out_dir, paste0(input_stem, "_plot_v2.png"))
svg_file <- file.path(out_dir, paste0(input_stem, "_plot_v2.svg"))

ggsave(pdf_file, p, width = fig_width, height = fig_height, units = "in", device = cairo_pdf)
ggsave(png_file, p, width = fig_width, height = fig_height, units = "in", dpi = 400, bg = "white")
ggsave(svg_file, p, width = fig_width, height = fig_height, units = "in")

print(p)
message("PDF: ", pdf_file)
message("PNG: ", png_file)
message("SVG: ", svg_file)
