#!/usr/bin/env Rscript

packages <- c("ggplot2", "dplyr", "readr", "tidyr", "patchwork")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(patchwork)

message("请选择包含 PT2906 和 U. ramanniana BUSCO 结果的 CSV 文件")
input_file <- file.choose()
out_dir <- dirname(input_file)
input_stem <- tools::file_path_sans_ext(basename(input_file))

raw <- read_csv(input_file, show_col_types = FALSE, trim_ws = TRUE)
required_cols <- c("sample", "dataset", "metric", "percent")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("数据缺少以下列：", paste(missing_cols, collapse = ", "))
}

clean_data <- raw %>%
  mutate(
    source_row = row_number(),
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
    sample %in% c("Rhizopus_sp._PT2906", "Umbelopsis_ramanniana"),
    dataset %in% c("Genome", "Longest-isoform protein"),
    metric %in% c("S", "D", "F", "M"),
    !is.na(percent)
  ) %>%
  arrange(source_row) %>%
  group_by(sample, dataset, metric) %>%
  mutate(run_index = row_number()) %>%
  ungroup() %>%
  mutate(
    assay = case_when(
      dataset == "Genome" ~ "Genome",
      sample == "Rhizopus_sp._PT2906" & dataset == "Longest-isoform protein" ~
        "Longest-isoform protein",
      sample == "Umbelopsis_ramanniana" &
        dataset == "Longest-isoform protein" & run_index == 1 ~
        "BRAKER3 longest-isoform protein",
      sample == "Umbelopsis_ramanniana" &
        dataset == "Longest-isoform protein" & run_index == 2 ~
        "BRAKER2 longest-isoform protein",
      TRUE ~ NA_character_
    ),
    category = recode(
      metric,
      S = "Complete (Single-copy)",
      D = "Complete (Duplicated)",
      F = "Fragmented",
      M = "Missing"
    )
  ) %>%
  filter(!is.na(assay))

expected_counts <- tibble(
  sample = c("Rhizopus_sp._PT2906", "Umbelopsis_ramanniana"),
  expected_rows = c(8L, 12L)
)
observed_counts <- clean_data %>% count(sample, name = "observed_rows")
count_check <- expected_counts %>%
  left_join(observed_counts, by = "sample") %>%
  mutate(observed_rows = replace_na(observed_rows, 0L))

if (any(count_check$expected_rows != count_check$observed_rows)) {
  print(count_check)
  stop(
    "输入行数不符合预期：PT2906应为8行；U. ramanniana应为12行。",
    "每个Genome或蛋白版本均须各有S、D、F、M四行。"
  )
}

combination_check <- clean_data %>%
  count(sample, assay, metric) %>%
  filter(n != 1)
if (nrow(combination_check) > 0) {
  print(combination_check)
  stop("部分 sample-assay-metric 组合不是唯一记录，请检查CSV行顺序和重复记录。")
}

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

category_levels <- c(
  "Missing",
  "Fragmented",
  "Complete (Duplicated)",
  "Complete (Single-copy)"
)

clean_data <- clean_data %>%
  mutate(category = factor(category, levels = category_levels))

bar_totals <- clean_data %>%
  group_by(sample, assay) %>%
  summarise(raw_total = sum(percent), .groups = "drop") %>%
  mutate(
    gap_to_100 = 100 - raw_total,
    exceeds_100 = raw_total > 100.05
  )

if (any(bar_totals$exceeds_100)) {
  warning("部分数据条S+D+F+M超过100%；脚本不会缩放原始数值。")
}

check_file <- file.path(out_dir, paste0(input_stem, "_panels_ab_v2_data_check.csv"))
write_csv(bar_totals, check_file)

make_panel <- function(data, assay_levels, assay_labels, panel_title) {
  panel_data <- data %>%
    mutate(assay = factor(assay, levels = assay_levels))

  label_data <- panel_data %>% filter(percent >= 3)

  ggplot(panel_data, aes(x = assay, y = percent, fill = category)) +
    geom_col(width = 0.70, colour = NA, position = "stack") +
    geom_text(
      data = label_data,
      aes(label = sprintf("%.1f", percent)),
      position = position_stack(vjust = 0.5),
      colour = "black",
      size = 3.15,
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
    scale_x_discrete(
      labels = assay_labels,
      drop = FALSE,
      expand = expansion(add = c(0.32, 0.32))
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = "BUSCO percentage (%)",
      fill = "BUSCO category"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = 10, base_family = "serif") +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold.italic",
        size = 11.5,
        margin = margin(b = 7)
      ),
      axis.title.x = element_text(face = "bold", size = 10.5),
      axis.text.x = element_text(size = 9, colour = "black"),
      axis.text.y = element_text(size = 9, colour = "black"),
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
}

panel_a_data <- clean_data %>% filter(sample == "Rhizopus_sp._PT2906")
panel_b_data <- clean_data %>% filter(sample == "Umbelopsis_ramanniana")

p_a <- make_panel(
  panel_a_data,
  assay_levels = c("Longest-isoform protein", "Genome"),
  assay_labels = c(
    "Genome" = "Genome (G)",
    "Longest-isoform protein" = "Protein (P)"
  ),
  panel_title = "Rhizopus sp. PT2906"
) +
  theme(plot.margin = margin(8, 12, 22, 8))

p_b <- make_panel(
  panel_b_data,
  assay_levels = c(
    "BRAKER2 longest-isoform protein",
    "BRAKER3 longest-isoform protein",
    "Genome"
  ),
  assay_labels = c(
    "Genome" = "Genome (G)",
    "BRAKER3 longest-isoform protein" = "BRAKER3 protein (B3)",
    "BRAKER2 longest-isoform protein" = "BRAKER2 protein (B2)"
  ),
  panel_title = "Umbelopsis ramanniana"
) +
  theme(plot.margin = margin(22, 12, 8, 8))

combined_plot <- (p_a / p_b) +
  plot_layout(
    heights = c(1.0, 1.22),
    guides = "collect"
  ) +
  plot_annotation(
    title = "BUSCO-guided resolution of quality-control outliers through reannotation and sample exclusion",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14,
        family = "serif",
        margin = margin(b = 8)
      ),
      plot.tag = element_text(
        face = "bold",
        size = 13,
        family = "serif"
      )
    )
  ) &
  theme(
    legend.position = "top",
    legend.justification = "center"
  )

pdf_file <- file.path(out_dir, paste0(input_stem, "_panels_ab_v2.pdf"))
png_file <- file.path(out_dir, paste0(input_stem, "_panels_ab_v2.png"))
svg_file <- file.path(out_dir, paste0(input_stem, "_panels_ab_v2.svg"))

ggsave(
  pdf_file, combined_plot,
  width = 9.2, height = 10.4,
  units = "in", device = cairo_pdf, limitsize = FALSE
)
ggsave(
  png_file, combined_plot,
  width = 9.2, height = 10.4,
  units = "in", dpi = 400,
  bg = "white", limitsize = FALSE
)
ggsave(
  svg_file, combined_plot,
  width = 9.2, height = 10.4,
  units = "in", limitsize = FALSE
)

print(combined_plot)
message("绘图完成。")
message("Panel a 位于上方：Rhizopus sp. PT2906。")
message("Panel b 位于下方：Umbelopsis ramanniana。")
message("两个panel之间保留空白；BUSCO四色与原v8脚本完全一致。")
message("PDF：", pdf_file)
message("PNG：", png_file)
message("SVG：", svg_file)
message("检查表：", check_file)
