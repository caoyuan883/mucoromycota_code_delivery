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
if (length(missing_cols) > 0) stop("数据缺少以下列：", paste(missing_cols, collapse = ", "))

sample_levels_original <- unique(trimws(as.character(raw$sample)))
sample_levels_original <- sample_levels_original[!is.na(sample_levels_original) & sample_levels_original != ""]

clean_data <- raw %>%
  mutate(
    sample = trimws(as.character(sample)),
    dataset = trimws(as.character(dataset)),
    metric = toupper(trimws(as.character(metric))),
    percent = as.numeric(gsub("%", "", as.character(percent))),
    dataset = case_when(
      dataset %in% c("G", "Genome") ~ "Genome",
      dataset %in% c("P", "Protein", "Longest-isoform protein") ~ "Longest-isoform protein",
      TRUE ~ dataset
    )
  ) %>%
  filter(
    sample != "", !is.na(sample), !is.na(percent),
    metric %in% c("S", "D", "F", "M"),
    dataset %in% c("Genome", "Longest-isoform protein")
  )

if (nrow(clean_data) == 0) stop("清理后没有可绘图的数据。")

expected_grid <- expand_grid(
  sample = sample_levels_original,
  dataset = c("Longest-isoform protein", "Genome"),
  metric = c("S", "D", "F", "M")
)

main_data <- expected_grid %>%
  left_join(
    clean_data %>%
      group_by(sample, dataset, metric) %>%
      summarise(percent = first(percent), .groups = "drop"),
    by = c("sample", "dataset", "metric")
  ) %>%
  mutate(
    percent = replace_na(percent, 0),
    dataset_label = ifelse(dataset == "Longest-isoform protein", "P", "G"),
    category = recode(
      metric,
      S = "Complete (Single-copy)",
      D = "Complete (Duplicated)",
      F = "Fragmented",
      M = "Missing"
    ),
    draw_order = recode(metric, S = 1L, D = 2L, F = 3L, M = 4L)
  )

bar_totals <- main_data %>%
  group_by(sample, dataset, dataset_label) %>%
  summarise(raw_total = sum(percent, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    gap = pmax(0, 100 - raw_total),
    exceeds_100 = raw_total > 100.05
  )

if (any(bar_totals$exceeds_100)) {
  warning("部分数据条S+D+F+M超过100%；未缩放原始值，请查看检查表。")
}

plot_data <- main_data %>%
  select(sample, dataset, dataset_label, category, draw_order, percent) %>%
  bind_rows(
    bar_totals %>%
      transmute(
        sample,
        dataset,
        dataset_label,
        category = "Rounding gap",
        draw_order = 5L,
        percent = gap
      )
  )

sample_levels_plot <- rev(sample_levels_original)
bar_levels <- unlist(
  lapply(sample_levels_plot, function(x) {
    c(paste0(x, "___Genome"), paste0(x, "___Longest-isoform protein"))
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
    bar = paste0(sample, "___", dataset),
    bar_index = match(bar, bar_levels)
  ) %>%
  arrange(bar_index, draw_order) %>%
  group_by(bar_index) %>%
  mutate(
    segment_start = cumsum(lag(percent, default = 0)),
    segment_end = segment_start + percent,
    segment_mid = segment_start + percent / 2
  ) %>%
  ungroup()

label_data <- plot_data %>%
  filter(category != "Rounding gap", percent >= 3) %>%
  mutate(percent_label = sprintf("%.1f", percent))

bar_frame_data <- tibble(
  bar_index = seq_along(bar_levels),
  xmin = 0,
  xmax = 100
)

check_file <- file.path(out_dir, paste0(input_stem, "_v6_data_check.csv"))
write_csv(bar_totals, check_file)

colors <- c(
  "Complete (Single-copy)" = "#2B8CBE",
  "Complete (Duplicated)" = "#084081",
  "Fragmented" = "#FDAE61",
  "Missing" = "#D7191C",
  "Rounding gap" = "#D9D9D9"
)

legend_order <- c(
  "Complete (Single-copy)",
  "Complete (Duplicated)",
  "Fragmented",
  "Missing"
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
bar_half_height <- 0.39

p <- ggplot() +
  geom_rect(
    data = bar_frame_data,
    aes(
      xmin = bar_index - bar_half_height,
      xmax = bar_index + bar_half_height,
      ymin = xmin,
      ymax = xmax
    ),
    inherit.aes = FALSE,
    fill = "white",
    colour = "black",
    linewidth = 0.38
  ) +
  geom_rect(
    data = plot_data %>% filter(percent > 0),
    aes(
      xmin = bar_index - bar_half_height,
      xmax = bar_index + bar_half_height,
      ymin = segment_start,
      ymax = segment_end,
      fill = category
    ),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 0.28
  ) +
  geom_rect(
    data = bar_frame_data,
    aes(
      xmin = bar_index - bar_half_height,
      xmax = bar_index + bar_half_height,
      ymin = xmin,
      ymax = xmax
    ),
    inherit.aes = FALSE,
    fill = NA,
    colour = "black",
    linewidth = 0.38
  ) +
  geom_text(
    data = label_data,
    aes(
      x = bar_index,
      y = segment_mid,
      label = percent_label
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.15,
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
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    breaks = seq_along(bar_levels),
    labels = axis_labels[bar_levels],
    limits = c(0.5, length(bar_levels) + 0.5),
    expand = expansion(mult = c(0, 0))
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
      colour = "grey85", linewidth = 0.35, linetype = "dashed"
    ),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(1.0, "cm"),
    plot.margin = margin(8, 12, 8, 8)
  )

pdf_file <- file.path(out_dir, paste0(input_stem, "_plot_v6.pdf"))
png_file <- file.path(out_dir, paste0(input_stem, "_plot_v6.png"))
svg_file <- file.path(out_dir, paste0(input_stem, "_plot_v6.svg"))

ggsave(
  pdf_file, p,
  width = fig_width, height = fig_height,
  units = "in", device = cairo_pdf, limitsize = FALSE
)

ggsave(
  png_file, p,
  width = fig_width, height = fig_height,
  units = "in", dpi = 400, bg = "white", limitsize = FALSE
)

ggsave(
  svg_file, p,
  width = fig_width, height = fig_height,
  units = "in", limitsize = FALSE
)

print(p)
message("绘图完成。")
message("每一条数据条都先画独立白色底框，再画色块，最后重画黑色外框。")
message("红色Missing即使为0，也保留整条100%白色底框，便于Adobe后续编辑。")
message("仅显示大于或等于3.0%的原始数值，全部数字为黑色。")
message("PDF：", pdf_file)
message("PNG：", png_file)
message("SVG：", svg_file)
message("检查表：", check_file)
