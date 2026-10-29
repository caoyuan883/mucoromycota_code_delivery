#!/usr/bin/env Rscript
# =============================================================================
# plot_busco_genome_vs_protein_by_subphylum_A4_v12.R
# BUSCO (mucoromycota_odb12): Genome (G) vs Longest-isoform protein (P)
# v12 相对 v11：
#   - 全部文字 Helvetica
#   - y 轴：属/种加词斜体；sp./var./菌株编号等正体；G / P 正体并与名称分开
#   - 柱内百分比粗体
#   - 删除全部图注（caption）
# 输入：longest_isoform_busco_genome_protein_mucoromycota_202_v3-2.xlsx
#       工作表 "1_简表_202株"
# =============================================================================

packages <- c("ggplot2", "dplyr", "tidyr", "readxl", "readr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(grid)

# ---------------- 参数 ----------------
sheet_name      <- "1_简表_202株"
name_column     <- "现用名"
subphylum_order <- c("Mucoromycotina", "Mortierellomycotina", "Glomeromycotina")
fig_letter      <- c(Mucoromycotina = "a", Mortierellomycotina = "b", Glomeromycotina = "c")

strains_per_page <- 41
bar_slot_in      <- 0.10
panel_width_mm   <- 115
max_fig_width_mm <- 180
bar_width        <- 0.82
label_min        <- 5

page_w_mm       <- 210
page_h_mm       <- 297
margin_top_mm   <- 15
margin_right_mm <- 15

show_title  <- TRUE
font_family <- "Helvetica"
fix_hyphen  <- TRUE   # R pdf 设备会把 "-" 画成减号；TRUE 时替换为 \u00AD 以输出真正的连字符

title_pt        <- 9
legend_title_pt <- 7.5
legend_text_pt  <- 7
axis_title_pt   <- 7.5
axis_x_text_pt  <- 7
axis_y_text_pt  <- 6.5
bar_label_pt    <- 5

fix_text <- function(x) if (fix_hyphen) gsub("-", "\u00AD", x, fixed = TRUE) else x

# ---------------- 读入 ----------------
message("请选择 longest_isoform_busco_genome_protein_mucoromycota_202_v3-2.xlsx")
input_file <- file.choose()
out_dir <- file.path(dirname(input_file), "BUSCO_G_vs_P_A4_v12")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
input_stem <- tools::file_path_sans_ext(basename(input_file))

raw <- read_excel(input_file, sheet = sheet_name)
required_cols <- c(name_column, "subphylum", "dataset", "metric", "percent")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) stop("数据缺少以下列：", paste(missing_cols, collapse = ", "))

# ---------------- 清理 ----------------
clean_data <- raw %>%
  transmute(
    name      = trimws(as.character(.data[[name_column]])),
    subphylum = trimws(as.character(subphylum)),
    dataset   = trimws(as.character(dataset)),
    metric    = toupper(trimws(as.character(metric))),
    percent   = as.numeric(gsub("%", "", as.character(percent)))
  ) %>%
  mutate(dataset = case_when(
    dataset %in% c("G", "Genome") ~ "Genome",
    dataset %in% c("P", "Protein", "Longest-isoform protein") ~ "Longest-isoform protein",
    TRUE ~ dataset
  )) %>%
  filter(!is.na(name), name != "", !is.na(percent),
         metric %in% c("S", "D", "F", "M"),
         dataset %in% c("Genome", "Longest-isoform protein"))

unknown_sub <- setdiff(unique(clean_data$subphylum), subphylum_order)
if (length(unknown_sub) > 0) stop("出现未定义的亚门：", paste(unknown_sub, collapse = ", "))
dup_check <- clean_data %>% distinct(name, subphylum) %>% count(name) %>% filter(n > 1)
if (nrow(dup_check) > 0) stop("以下名称对应多个亚门：", paste(dup_check$name, collapse = ", "))

sample_info <- clean_data %>%
  distinct(name, subphylum) %>%
  mutate(subphylum = factor(subphylum, levels = subphylum_order)) %>%
  arrange(subphylum, tolower(name))

plot_data <- sample_info %>%
  tidyr::crossing(dataset = c("Genome", "Longest-isoform protein"),
                  metric  = c("S", "D", "F", "M")) %>%
  left_join(clean_data %>% group_by(name, dataset, metric) %>%
              summarise(percent = first(percent), .groups = "drop"),
            by = c("name", "dataset", "metric")) %>%
  mutate(
    percent  = replace_na(percent, 0),
    category = recode(metric, S = "Complete (Single-copy)", D = "Complete (Duplicated)",
                      F = "Fragmented", M = "Missing"),
    category = factor(category, levels = c("Missing", "Fragmented",
                                           "Complete (Duplicated)", "Complete (Single-copy)"))
  )

bar_totals <- plot_data %>%
  group_by(subphylum, name, dataset) %>%
  summarise(raw_total = sum(percent), .groups = "drop") %>%
  mutate(gap_to_100 = round(100 - raw_total, 2), exceeds_100 = raw_total > 100.2)
if (any(bar_totals$exceeds_100)) warning("部分数据条 S+D+F+M 超过 100.2%，请查看检查表。")
write_csv(bar_totals, file.path(out_dir, paste0(input_stem, "_v12_data_check.csv")))

# ---------------- 分页 ----------------
split_even <- function(n, k) {
  sizes <- rep(n %/% k, k) + c(rep(1, n %% k), rep(0, k - n %% k))
  rep(seq_len(k), sizes)
}

sample_info <- sample_info %>%
  group_by(subphylum) %>%
  mutate(
    n_sub   = n(),
    n_parts = ceiling(n_sub / strains_per_page),
    part    = split_even(n_sub[1], n_parts[1]),
    idx     = row_number()
  ) %>%
  ungroup()

message("分页方案：")
print(sample_info %>% count(subphylum, part, name = "strains"))

# ---------------- 颜色 ----------------
colors <- c("Complete (Single-copy)" = "#2B8CBE", "Complete (Duplicated)" = "#084081",
            "Fragmented" = "#FDAE61", "Missing" = "#D7191C")
legend_order  <- names(colors)
legend_labels <- fix_text(legend_order)

# ---------------- y 轴标签（plotmath）----------------
# 规则：第 1 个词（属名）斜体；全小写词（种加词、变种加词）斜体；
#       含 "." 或大写/数字的词（sp., var., f., UoMAU379 等）正体；G/P 正体
make_axis_expr <- function(name, ds_short) {
  dn  <- gsub("_", " ", name, fixed = TRUE)
  tok <- strsplit(trimws(dn), "\\s+")[[1]]
  is_italic <- seq_along(tok) == 1 | grepl("^[a-z][a-z-]*$", tok)
  tok <- fix_text(tok)
  q   <- sprintf('"%s"', gsub('(["\\\\])', '\\\\\\1', tok))
  parts <- ifelse(is_italic, paste0("italic(", q, ")"), q)
  paste(c(parts, sprintf('"%s"', ds_short)), collapse = "~")
}

# ---------------- 作图 ----------------
make_plot <- function(info, title_text) {
  names_plot <- rev(info$name)
  bar_levels <- unlist(lapply(names_plot, function(x)
    c(paste0(x, "___Longest-isoform protein"), paste0(x, "___Genome"))), use.names = FALSE)
  label_str <- unlist(lapply(names_plot, function(x)
    c(make_axis_expr(x, "P"), make_axis_expr(x, "G"))), use.names = FALSE)
  axis_labels <- parse(text = label_str)

  d <- plot_data %>%
    filter(name %in% info$name) %>%
    mutate(bar = factor(paste0(name, "___", dataset), levels = bar_levels))

  lab <- d %>%
    filter(percent >= label_min) %>%
    group_by(bar) %>% arrange(category, .by_group = TRUE) %>%
    mutate(pos = cumsum(percent) - percent / 2, txt = sprintf("%.1f", percent)) %>%
    ungroup()

  p <- ggplot(d, aes(x = percent, y = bar, fill = category)) +
    geom_col(width = bar_width, colour = NA, position = "stack", orientation = "y") +
    geom_text(data = lab, aes(x = pos, y = bar, label = txt), inherit.aes = FALSE,
              colour = "black", size = bar_label_pt / .pt, fontface = "bold",
              family = font_family) +
    scale_fill_manual(values = colors, breaks = legend_order, labels = legend_labels,
                      drop = FALSE) +
    scale_x_continuous(breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(breaks = bar_levels, labels = axis_labels, drop = FALSE,
                     expand = expansion(add = c(0.5, 0.5))) +
    coord_cartesian(xlim = c(0, 100), clip = "on") +
    labs(title = if (show_title) title_text else NULL,
         x = "BUSCO percentage (%)", y = NULL, fill = "BUSCO category") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = 7, base_family = font_family) +
    theme(
      text               = element_text(family = font_family),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = title_pt),
      axis.title.x       = element_text(face = "bold", size = axis_title_pt),
      axis.text.x        = element_text(size = axis_x_text_pt, colour = "black"),
      axis.text.y        = element_text(size = axis_y_text_pt, face = "plain", colour = "black"),
      axis.ticks         = element_line(linewidth = 0.25),
      panel.border       = element_rect(linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey85", linewidth = 0.25, linetype = "dashed"),
      legend.position    = "top",
      legend.title       = element_text(face = "bold", size = legend_title_pt),
      legend.text        = element_text(size = legend_text_pt),
      legend.key.size    = unit(3, "mm"),
      legend.key.width   = unit(6, "mm"),
      legend.margin      = margin(0, 0, 0, 0),
      plot.margin        = margin(2, 4, 2, 2)
    )
  p
}

fix_panel <- function(p, n_bars) {
  g  <- ggplotGrob(p)
  pr <- g$layout[grepl("^panel", g$layout$name), ]
  g$heights[unique(pr$t)] <- unit(n_bars * bar_slot_in, "in")
  g$widths[unique(pr$l)]  <- unit(panel_width_mm, "mm")
  list(g = g,
       w = convertWidth(sum(g$widths), "in", valueOnly = TRUE),
       h = convertHeight(sum(g$heights), "in", valueOnly = TRUE))
}

mm2in  <- function(x) x / 25.4
page_w <- mm2in(page_w_mm); page_h <- mm2in(page_h_mm)

parts <- sample_info %>% distinct(subphylum, part, n_parts, n_sub) %>% arrange(subphylum, part)
size_log <- list()

for (r in seq_len(nrow(parts))) {
  sp <- as.character(parts$subphylum[r]); k <- parts$part[r]
  n_parts <- parts$n_parts[r]; n_sub <- parts$n_sub[r]
  info <- sample_info %>% filter(subphylum == sp, part == k)
  i_from <- min(info$idx); i_to <- max(info$idx)

  title_text <- if (n_parts > 1) {
    sprintf("%s \u00B7 part %d/%d \u00B7 strains %d\u2013%d of %d", sp, k, n_parts, i_from, i_to, n_sub)
  } else {
    sprintf("%s \u00B7 n = %d", sp, n_sub)
  }

  tag <- if (n_parts > 1) sprintf("%s%d", fig_letter[[sp]], k) else fig_letter[[sp]]
  f <- file.path(out_dir, sprintf("Fig_S%s_%s%s_v12.pdf", tag, sp,
                                  if (n_parts > 1) sprintf("_part%d", k) else ""))
  if (file.exists(f)) file.remove(f)

  grDevices::pdf(f, width = page_w, height = page_h, family = font_family, useDingbats = FALSE)
  grid.newpage()
  fp <- fix_panel(make_plot(info, title_text), nrow(info) * 2)

  if (fp$w > mm2in(max_fig_width_mm)) {
    warning(sprintf("%s：图宽 %.1f mm 超过上限 %d mm，请减小 panel_width_mm 或 axis_y_text_pt",
                    basename(f), fp$w * 25.4, max_fig_width_mm))
  }

  x_right <- page_w - mm2in(margin_right_mm)
  y_top   <- page_h - mm2in(margin_top_mm)
  pushViewport(viewport(x = unit(x_right, "in"), y = unit(y_top, "in"),
                        width = unit(fp$w, "in"), height = unit(fp$h, "in"),
                        just = c("right", "top")))
  grid.draw(fp$g)
  popViewport()
  grDevices::dev.off()

  if (!(file.exists(f) && file.info(f)$size > 0)) stop("PDF 未生成：", f)
  size_log[[r]] <- data.frame(file = basename(f), strains = nrow(info),
                              fig_width_mm  = round(fp$w * 25.4, 1),
                              fig_height_mm = round(fp$h * 25.4, 1))
  message("已输出：", f)
}

size_tab <- do.call(rbind, size_log)
print(size_tab)
write_csv(size_tab, file.path(out_dir, "figure_sizes_v12.csv"))
message("完成。输出目录：", out_dir)
