#!/usr/bin/env Rscript

packages <- c("ggplot2", "dplyr", "readr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(readr)

message("请选择 longest_isoform_busco_genome_protein_*_part*.csv")
input_file <- file.choose()
out_dir <- dirname(input_file)
input_stem <- tools::file_path_sans_ext(basename(input_file))

raw <- read_csv(input_file, show_col_types = FALSE, trim_ws = TRUE)
required_cols <- c("sample", "dataset", "dataset_label", "metric", "percent")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("数据缺少以下列：", paste(missing_cols, collapse = ", "))
}

raw$input_row <- seq_len(nrow(raw))
sample_levels_original <- unique(trimws(as.character(raw$sample)))
sample_levels_original <- sample_levels_original[
  !is.na(sample_levels_original) & sample_levels_original != ""
]

plot_data <- raw %>%
  mutate(
    sample = trimws(as.character(sample)),
    dataset = trimws(as.character(dataset)),
    metric = toupper(trimws(as.character(metric))),
    percent = as.numeric(gsub("%", "", as.character(percent))),
    dataset = case_when(
      dataset %in% c("G", "Genome") ~ "Genome",
      dataset %in% c("P", "Protein", "Longest-isoform protein") ~ "Longest-isoform protein",
      TRUE ~ dataset
    ),
    dataset = factor(dataset, levels = c("Genome", "Longest-isoform protein")),
    dataset_label = recode(
      as.character(dataset),
      "Genome" = "G",
      "Longest-isoform protein" = "P"
    ),
    metric = factor(metric, levels = c("S", "D", "F", "M")),
    category = recode(
      as.character(metric),
      S = "Complete (Single-copy)",
      D = "Complete (Duplicated)",
      F = "Fragmented",
      M = "Missing"
    ),
    category = factor(
      category,
      levels = c(
        "Missing",
        "Fragmented",
        "Complete (Duplicated)",
        "Complete (Single-copy)"
      )
    )
  ) %>%
  filter(sample != "", !is.na(sample), !is.na(percent), !is.na(metric), !is.na(dataset))

if (nrow(plot_data) == 0) stop("清理后没有可绘图的数据。")

check <- plot_data %>%
  group_by(sample, dataset) %>%
  summarise(
    n_category = n_distinct(metric),
    total = sum(percent, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(check, file.path(out_dir, paste0(input_stem, "_data_check.csv")))

bad_check <- check %>% filter(n_category != 4 | total < 99.5 | total > 100.5)
if (nrow(bad_check) > 0) {
  warning("部分样本缺少类别或S+D+F+M不接近100%；请查看_data_check.csv。")
}

sample_levels_plot <- rev(sample_levels_original)

bar_levels <- unlist(
  lapply(sample_levels_plot, function(x) {
    c(
      paste0(x, "___Longest-isoform protein"),
      paste0(x, "___Genome")
    )
  }),
  use.names = FALSE
)

species_label_values <- unlist(
  lapply(sample_levels_plot, function(x) c(x, "")),
  use.names = FALSE
)
species_labels <- setNames(species_label_values, bar_levels)

plot_data <- plot_data %>%
  mutate(
    sample = factor(sample, levels = sample_levels_plot),
    bar = paste0(as.character(sample), "___", as.character(dataset)),
    bar = factor(bar, levels = bar_levels)
  )

colors <- c(
  "Complete (Single-copy)" = "#2B8CBE",
  "Complete (Duplicated)" = "#084081",
  "Fragmented" = "#FDAE61",
  "Missing" = "#D7191C"
)

legend_order <- c(
  "Complete (Single-copy)",
  "Complete (Duplicated)",
  "Fragmented",
  "Missing"
)

label_data <- plot_data %>%
  distinct(bar, dataset_label) %>%
  mutate(y = 100.6)

lineage_label <- if ("lineage" %in% names(plot_data)) {
  unique(as.character(plot_data$lineage))[1]
} else {
  "BUSCO"
}

lineage_title <- ifelse(
  tolower(lineage_label) == "mucoromycota",
  "Mucoromycota BUSCO",
  ifelse(tolower(lineage_label) == "fungi", "Fungi BUSCO", "BUSCO")
)

n_samples <- length(sample_levels_original)
fig_width <- 13
fig_height <- max(8, n_samples * 0.27 + 2.5)

p <- ggplot(plot_data, aes(x = bar, y = percent, fill = category)) +
  geom_col(width = 0.82, colour = NA, position = "stack") +
  geom_text(
    data = label_data,
    aes(x = bar, y = y, label = dataset_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.5,
    fontface = "bold",
    family = "serif"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = colors,
    breaks = legend_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 102),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    labels = species_labels,
    drop = FALSE,
    expand = expansion(add = c(0.25, 0.25))
  ) +
  labs(
    title = paste0(
      "Comparison of ", lineage_title,
      " Scores for Genome (G) and Longest-isoform Protein (P)"
    ),
    x = NULL,
    y = "BUSCO percentage (%)",
    fill = "BUSCO category"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_bw(base_size = 10, base_family = "serif") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12.5),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 7.5, face = "italic", colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      colour = "grey85",
      linewidth = 0.35,
      linetype = "dashed"
    ),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(1.1, "cm"),
    plot.margin = margin(8, 25, 8, 8)
  )

pdf_file <- file.path(out_dir, paste0(input_stem, "_plot_v2.pdf"))
png_file <- file.path(out_dir, paste0(input_stem, "_plot_v2.png"))
svg_file <- file.path(out_dir, paste0(input_stem, "_plot_v2.svg"))

ggsave(
  pdf_file, p,
  width = fig_width, height = fig_height,
  units = "in", device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  png_file, p,
  width = fig_width, height = fig_height,
  units = "in", dpi = 400,
  bg = "white", limitsize = FALSE
)

ggsave(
  svg_file, p,
  width = fig_width, height = fig_height,
  units = "in", limitsize = FALSE
)

print(p)
message("绘图完成。图中从上到下与CSV中的物种顺序一致。")
message("每株从上到下：G、P。")
message("每条柱从左到右：Single-copy、Duplicated、Fragmented、Missing。")
message("PDF：", pdf_file)
message("PNG：", png_file)
message("SVG：", svg_file)
