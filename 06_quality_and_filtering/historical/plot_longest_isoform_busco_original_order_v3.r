#!/usr/bin/env Rscript

packages <- c("ggplot2", "dplyr", "readr", "tidyr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

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
    dataset_label = case_when(
      dataset == "Genome" ~ "G",
      dataset == "Longest-isoform protein" ~ "P",
      TRUE ~ as.character(dataset_label)
    ),
    metric = factor(metric, levels = c("S", "D", "F", "M")),
    category = recode(
      as.character(metric),
      S = "Complete (Single-copy)",
      D = "Complete (Duplicated)",
      F = "Fragmented",
      M = "Missing"
    )
  ) %>%
  filter(
    sample != "",
    !is.na(sample),
    !is.na(percent),
    !is.na(metric),
    dataset %in% c("Genome", "Longest-isoform protein")
  )

if (nrow(plot_data) == 0) stop("清理后没有可绘图的数据。")

expected_categories <- c("S", "D", "F", "M")
expected_grid <- expand_grid(
  sample = sample_levels_original,
  dataset = c("Longest-isoform protein", "Genome"),
  metric_chr = expected_categories
)

plot_data <- expected_grid %>%
  left_join(
    plot_data %>%
      transmute(
        sample,
        dataset,
        metric_chr = as.character(metric),
        percent
      ),
    by = c("sample", "dataset", "metric_chr")
  ) %>%
  mutate(
    percent = replace_na(percent, 0),
    dataset_label = ifelse(dataset == "Longest-isoform protein", "P", "G"),
    category = recode(
      metric_chr,
      S = "Complete (Single-copy)",
      D = "Complete (Duplicated)",
      F = "Fragmented",
      M = "Missing"
    )
  )

bar_totals <- plot_data %>%
  group_by(sample, dataset, dataset_label) %>%
  summarise(raw_total = sum(percent, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    unclassified = pmax(0, 100 - raw_total),
    exceeds_100 = raw_total > 100.05
  )

if (any(bar_totals$exceeds_100)) {
  warning("部分数据条S+D+F+M超过100%；未缩放原始值，请查看检查表。")
}

plot_data <- plot_data %>%
  select(sample, dataset, dataset_label, category, percent) %>%
  bind_rows(
    bar_totals %>%
      transmute(
        sample,
        dataset,
        dataset_label,
        category = "Unclassified / rounding gap",
        percent = unclassified
      )
  )

stack_levels <- c(
  "Unclassified / rounding gap",
  "Missing",
  "Fragmented",
  "Complete (Duplicated)",
  "Complete (Single-copy)"
)

plot_data <- plot_data %>%
  mutate(category = factor(category, levels = stack_levels))

sample_levels_plot <- rev(sample_levels_original)
bar_levels <- unlist(
  lapply(sample_levels_plot, function(x) {
    c(
      paste0(x, "___Genome"),
      paste0(x, "___Longest-isoform protein")
    )
  }),
  use.names = FALSE
)

axis_labels <- setNames(
  unlist(
    lapply(sample_levels_plot, function(x) {
      c(paste0(x, "  G"), paste0(x, "  P"))
    }),
    use.names = FALSE
  ),
  bar_levels
)

plot_data <- plot_data %>%
  mutate(
    sample = factor(sample, levels = sample_levels_plot),
    bar = paste0(as.character(sample), "___", dataset),
    bar = factor(bar, levels = bar_levels)
  )

label_data <- plot_data %>%
  filter(percent > 0) %>%
  group_by(bar) %>%
  arrange(category, .by_group = TRUE) %>%
  mutate(
    label_position = cumsum(percent) - percent / 2,
    percent_label = ifelse(
      as.character(category) == "Unclassified / rounding gap",
      sprintf("%.1f", percent),
      sprintf("%.1f", percent)
    ),
    text_colour = case_when(
      as.character(category) %in% c(
        "Complete (Single-copy)",
        "Complete (Duplicated)",
        "Missing"
      ) ~ "white",
      TRUE ~ "black"
    ),
    show_label = percent >= 1.5
  ) %>%
  ungroup() %>%
  filter(show_label)

small_label_data <- plot_data %>%
  filter(
    percent > 0,
    percent < 1.5,
    as.character(category) != "Unclassified / rounding gap"
  ) %>%
  group_by(bar) %>%
  arrange(category, .by_group = TRUE) %>%
  mutate(
    cumulative_end = cumsum(percent),
    percent_label = sprintf("%.1f", percent)
  ) %>%
  ungroup()

check_file <- file.path(out_dir, paste0(input_stem, "_v3_data_check.csv"))
write_csv(bar_totals, check_file)

colors <- c(
  "Complete (Single-copy)" = "#2B8CBE",
  "Complete (Duplicated)" = "#084081",
  "Fragmented" = "#FDAE61",
  "Missing" = "#D7191C",
  "Unclassified / rounding gap" = "#D9D9D9"
)

legend_order <- c(
  "Complete (Single-copy)",
  "Complete (Duplicated)",
  "Fragmented",
  "Missing",
  "Unclassified / rounding gap"
)

lineage_label <- if ("lineage" %in% names(raw)) {
  unique(trimws(as.character(raw$lineage)))[1]
} else {
  "BUSCO"
}
lineage_title <- ifelse(
  tolower(lineage_label) == "mucoromycota",
  "Mucoromycota BUSCO",
  ifelse(tolower(lineage_label) == "fungi", "Fungi BUSCO", "BUSCO")
)

n_samples <- length(sample_levels_original)
fig_width <- 16
fig_height <- max(9, n_samples * 0.38 + 2.8)

p <- ggplot(plot_data, aes(x = bar, y = percent, fill = category)) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.12, position = "stack") +
  geom_text(
    data = label_data,
    aes(x = bar, y = label_position, label = percent_label, colour = text_colour),
    inherit.aes = FALSE,
    size = 2.15,
    fontface = "bold",
    family = "serif"
  ) +
  scale_colour_identity() +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = colors,
    breaks = legend_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    labels = axis_labels,
    drop = FALSE,
    expand = expansion(add = c(0.22, 0.22))
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
    axis.text.y = element_text(size = 6.8, face = "italic", colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      colour = "grey85",
      linewidth = 0.35,
      linetype = "dashed"
    ),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(1.0, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  )

pdf_file <- file.path(out_dir, paste0(input_stem, "_plot_v3.pdf"))
png_file <- file.path(out_dir, paste0(input_stem, "_plot_v3.png"))
svg_file <- file.path(out_dir, paste0(input_stem, "_plot_v3.svg"))

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
message("绘图完成。")
message("每株从上到下：P、G；两条均显示物种名和P/G标记。")
message("不足100%的部分已用灰色补齐。")
message("大于或等于1.5%的区段显示原始百分比；更小区段保留在图中但不标字，以避免重叠。")
message("PDF：", pdf_file)
message("PNG：", png_file)
message("SVG：", svg_file)
message("检查表：", check_file)
