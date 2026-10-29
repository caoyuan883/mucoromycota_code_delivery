## =============================================================================
## Fig. 2  Annotation-protein species tree (technical validation) — v3
##  a: P 树(v7 画法) + 与基因组树的差异标记 + BRAKER 圈 + InterProScan / eggNOG-mapper 圈
##     + 外圈三段柱
##  b: 一张小提琴图，三个亚门三列（两种工具合计的注释比例）
## 与 v2 相同的输入文件（同一文件夹）：
##   P_occ80_g90.treefile  Fig2_tip_data.tsv  Fig2_conflict_splits.tsv  Fig2_conflict_regions.tsv
## =============================================================================

## ---- 0. settings ------------------------------------------------------------
label_size       <- 1.5
shade_alpha      <- 0.14
open_angle       <- 12
label_band_scale <- 1.00
ring_w           <- 0.035
ring_gap         <- 0.008
bar_h            <- 0.20
mark_conflict_branches <- TRUE
draw_region_arcs       <- TRUE
comparison_note  <- "Shared splits with genome tree: 197/204; nRF = 0.034 (outgroups retained)"
cov_limits       <- c(40, 100)
violin_var       <- "combined_pct"          # 也可改为 "interproscan_pct" 或 "eggnog_pct"
violin_ylab      <- "Annotated proteins (%)"
panel_b_position <- "inset"                 # "inset" = 放在树图右下角；"bottom" = 放在树图下方
fig_width_in     <- 14
fig_height_in    <- 13
out_prefix       <- "Fig2_protein_tree_annotation"

## ---- 1. packages ------------------------------------------------------------
pkgs_cran <- c("ape","ggplot2","dplyr","readr","stringr","ggnewscale","patchwork","scales","rstudioapi")
pkgs_bioc <- c("ggtree","treeio")
need <- setdiff(c(pkgs_cran, pkgs_bioc), rownames(installed.packages()))
if (length(need)) {
  if (length(intersect(need, pkgs_cran))) install.packages(intersect(need, pkgs_cran))
  if (length(intersect(need, pkgs_bioc))) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(intersect(need, pkgs_bioc), update = FALSE, ask = FALSE)
  }
}
suppressPackageStartupMessages({
  library(ape); library(treeio); library(ggtree); library(ggplot2); library(ggnewscale)
  library(dplyr); library(readr); library(stringr); library(patchwork); library(scales)
})

## ---- 2. files ---------------------------------------------------------------
work_dir <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) "")
has_all <- function(d) nzchar(d) && dir.exists(d) &&
  length(list.files(d, pattern = "^P_.*treefile$")) && file.exists(file.path(d, "Fig2_tip_data.tsv"))
if (!has_all(work_dir)) {
  message("请在弹窗中选中 Fig2_tip_data.tsv（其余文件须在同一文件夹）")
  work_dir <- dirname(file.choose())
}
if (!has_all(work_dir)) stop("文件夹里缺少 P_*.treefile 或 Fig2_tip_data.tsv: ", work_dir)
setwd(work_dir)
p_file     <- list.files(work_dir, pattern = "^P_.*treefile$", full.names = TRUE)[1]
meta_file  <- file.path(work_dir, "Fig2_tip_data.tsv")
split_file <- file.path(work_dir, "Fig2_conflict_splits.tsv")
reg_file   <- file.path(work_dir, "Fig2_conflict_regions.tsv")

## ---- 3. read, root ------------------------------------------------------------
tr <- read.tree(p_file)
tr$tip.label <- gsub("^'|'$", "", trimws(tr$tip.label))
meta <- read_tsv(meta_file, show_col_types = FALSE, na = c("", "NA"))
bad1 <- setdiff(tr$tip.label, meta$label); bad2 <- setdiff(meta$label, tr$tip.label)
if (length(bad1) || length(bad2)) stop("名字对不上。树有表无: ", paste(bad1, collapse = ", "),
                                       "；表有树无: ", paste(bad2, collapse = ", "))
og <- meta$label[meta$role == "outgroup"]
if (!is.monophyletic(tr, og)) stop("外类群不成一支，不能定根")
tr <- root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
tr <- ladderize(tr, right = FALSE)

nl <- if (is.null(tr$node.label)) rep("", tr$Nnode) else tr$node.label
sp <- str_split_fixed(nl, "/", 2)
sup <- tibble(node = Ntip(tr) + seq_along(nl),
              shalrt = suppressWarnings(as.numeric(sp[, 1])),
              ufboot = suppressWarnings(as.numeric(sp[, 2]))) %>%
  mutate(weak = !is.na(shalrt) & !is.na(ufboot) & (shalrt < 80 | ufboot < 95))
weak_nodes <- sup$node[sup$weak]

conf_nodes <- integer(0)
if (file.exists(split_file)) {
  cs_tab <- read_tsv(split_file, show_col_types = FALSE)
  for (i in seq_len(nrow(cs_tab))) {
    tp <- strsplit(cs_tab$split_tips[i], ";")[[1]]
    if (!all(tp %in% tr$tip.label) || !is.monophyletic(tr, tp)) { warning("无法定位: ", cs_tab$split_tips[i]); next }
    conf_nodes <- c(conf_nodes, getMRCA(tr, tp))
  }
}

tr_plot <- compute.brlen(tr, method = "Grafen", power = 0.7)
tr_plot$node.label <- tr$node.label

## ---- 4. names ----------------------------------------------------------------
fmt_name <- function(x) {
  x <- gsub("\\s+", " ", trimws(x))
  x <- gsub("\\b(sp|cf|aff)\\.(?=\\S)", "\\1. ", x, perl = TRUE)
  esc <- function(s) gsub('"', '\\\\"', s)
  m <- regmatches(x, regexec("^(\\S+)\\s+(sp\\.|cf\\.|aff\\.)\\s*(.*)$", x))[[1]]
  if (length(m)) {
    sprintf('bolditalic("%s")~bold("%s")', esc(m[2]), esc(trimws(paste(m[3], m[4]))))
  } else {
    sprintf('bolditalic("%s")', esc(x))
  }
}
meta2 <- meta %>%
  mutate(display_name = ifelse(is.na(display_name), gsub("_", " ", label), display_name),
         group = ifelse(role == "outgroup", "Outgroup", group),
         braker_mode = ifelse(role == "outgroup", "Reference annotation", braker_mode),
         name_expr = vapply(display_name, fmt_name, character(1)))

## ---- 5. colours (v7) -----------------------------------------------------------
grp_levels <- c("Glomeromycotina", "Mortierellomycotina", "Mucoromycotina", "Outgroup")
grp_cols  <- c(Glomeromycotina = "#4DAF4A", Mortierellomycotina = "#984EA3", Mucoromycotina = "#377EB8", Outgroup = "grey55")
mode_cols <- c(BRAKER2 = "#E69F00", BRAKER3 = "#D55E00", "Reference annotation" = "grey70")
cov_cols  <- c("#FFFFE5", "#C7E9B4", "#41B6C4", "#225EA8", "#0C2C84")
seg_cols  <- c("Both tools" = "#4D4D4D", "eggNOG-mapper only" = "#56B4E9", "InterProScan only" = "#CC79A7")
conf_col  <- "#D7191C"

hl <- bind_rows(lapply(grp_levels, function(g) {
  tp <- meta2$label[meta2$group == g]
  if (!length(tp)) return(NULL)
  node <- if (length(tp) == 1) which(tr_plot$tip.label == tp) else getMRCA(tr_plot, tp)
  data.frame(lineage = g, node = node)
}))

## ---- 6. geometry --------------------------------------------------------------
p0 <- ggtree(tr_plot, layout = "fan", open.angle = open_angle, size = 0.3)
xmax <- max(p0$data$x[p0$data$isTip], na.rm = TRUE)
outer_pad <- ring_gap * 4 + ring_w * 3 + bar_h + 0.25
max_chars <- max(nchar(meta2$display_name))
label_in  <- max_chars * label_size * 0.62 / 25.4 * label_band_scale
panel_r_in <- 0.46 * min(fig_width_in - 2.8, fig_height_in)
ratio <- min(0.85, label_in * (1 + outer_pad) / panel_r_in)
R_in <- xmax; R_out <- R_in / (1 - ratio); R_text <- R_out - 0.008 * R_out
R_mode0 <- R_out  + ring_gap * xmax;        R_mode1 <- R_mode0 + ring_w * xmax
R_ips0  <- R_mode1 + ring_gap * xmax;       R_ips1  <- R_ips0  + ring_w * xmax
R_egg0  <- R_ips1  + ring_gap * xmax;       R_egg1  <- R_egg0  + ring_w * xmax
R_bar0  <- R_egg1  + 1.5 * ring_gap * xmax; H_bar   <- bar_h * xmax
R_arc   <- R_bar0  + H_bar + 0.03 * xmax

tips <- p0$data %>% filter(isTip) %>% select(label, y, angle) %>%
  left_join(meta2, by = "label") %>%
  mutate(flip = angle > 90 & angle < 270, tangle = ifelse(flip, angle + 180, angle), hj = ifelse(flip, 0, 1))
ing <- tips %>% filter(role == "ingroup")

## ---- 7. tree + lineage sector -------------------------------------------------
p <- p0 +
  geom_hilight(data = hl, mapping = aes(node = node, fill = lineage), alpha = shade_alpha) +
  geom_rect(data = tips, inherit.aes = FALSE,
            aes(xmin = R_in, xmax = R_out, ymin = y - 0.5, ymax = y + 0.5, fill = group),
            alpha = shade_alpha, colour = NA) +
  scale_fill_manual(values = grp_cols, breaks = grp_levels, name = "Lineage")
if (mark_conflict_branches && length(conf_nodes)) {
  pd <- p0$data
  cs <- pd %>% filter(node %in% conf_nodes) %>% mutate(xp = pd$x[match(parent, pd$node)])
  p <- p + geom_segment(data = cs, inherit.aes = FALSE,
                        aes(x = xp, xend = x, y = y, yend = y), colour = conf_col, linewidth = 0.8)
}
p <- p +
  geom_point2(aes(subset = node %in% weak_nodes), shape = 21, size = 1.4, fill = "white", colour = "black", stroke = 0.4) +
  geom_text(data = tips, inherit.aes = FALSE, parse = TRUE,
            aes(x = R_text, y = y, label = name_expr, angle = tangle, hjust = hj),
            size = label_size, colour = "grey10")

## ---- 8. rings + outer bars -----------------------------------------------------
p <- p + new_scale_fill() +
  geom_rect(data = tips, inherit.aes = FALSE,
            aes(xmin = R_mode0, xmax = R_mode1, ymin = y - 0.5, ymax = y + 0.5, fill = braker_mode), colour = NA) +
  scale_fill_manual(values = mode_cols, name = "Annotation")
p <- p + new_scale_fill() +
  geom_rect(data = ing, inherit.aes = FALSE,
            aes(xmin = R_ips0, xmax = R_ips1, ymin = y - 0.5, ymax = y + 0.5, fill = interproscan_pct), colour = NA) +
  geom_rect(data = ing, inherit.aes = FALSE,
            aes(xmin = R_egg0, xmax = R_egg1, ymin = y - 0.5, ymax = y + 0.5, fill = eggnog_pct), colour = NA) +
  scale_fill_gradientn(colours = cov_cols, limits = cov_limits, oob = squish,
                       name = "Annotated proteins (%)\ninner ring: InterProScan\nouter ring: eggNOG-mapper")
bars <- ing %>%
  transmute(y, s0 = R_bar0,
            s1 = s0 + H_bar * both_tools_pct / 100,
            s2 = s1 + H_bar * eggnog_only_pct / 100,
            s3 = s2 + H_bar * interproscan_only_pct / 100)
bars_long <- bind_rows(
  bars %>% transmute(y, xmin = s0, xmax = s1, seg = "Both tools"),
  bars %>% transmute(y, xmin = s1, xmax = s2, seg = "eggNOG-mapper only"),
  bars %>% transmute(y, xmin = s2, xmax = s3, seg = "InterProScan only")) %>%
  mutate(seg = factor(seg, levels = names(seg_cols)))
yr <- range(ing$y)
ref <- expand.grid(y = seq(yr[1] - 0.5, yr[2] + 0.5, length.out = 300), lvl = c(0.5, 1)) %>%
  mutate(x = R_bar0 + H_bar * lvl)
p <- p +
  geom_path(data = ref, inherit.aes = FALSE, aes(x = x, y = y, group = lvl),
            colour = "grey75", linewidth = 0.2, linetype = "dashed") +
  new_scale_fill() +
  geom_rect(data = bars_long, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = y - 0.4, ymax = y + 0.4, fill = seg), colour = NA) +
  scale_fill_manual(values = seg_cols, name = "Outer bars:\nannotated proteins (%)\ndashed lines: 50 %, 100 %")

## ---- 9. conflict region arcs -----------------------------------------------------
if (draw_region_arcs && file.exists(reg_file)) {
  rg <- read_tsv(reg_file, show_col_types = FALSE)
  arcs <- bind_rows(lapply(seq_len(nrow(rg)), function(i) {
    tp <- strsplit(rg$tips[i], ";")[[1]]
    sel <- tips$label %in% tp
    data.frame(region = rg$region[i], ymin = min(tips$y[sel]) - 0.5, ymax = max(tips$y[sel]) + 0.5,
               ymid = mean(range(tips$y[sel])), ang = mean(tips$angle[sel]))
  })) %>% mutate(flip = ang > 90 & ang < 270, tang = ifelse(flip, ang + 180, ang), hj = ifelse(flip, 1, 0))
  p <- p +
    geom_rect(data = arcs, inherit.aes = FALSE,
              aes(xmin = R_arc, xmax = R_arc + 0.006 * xmax, ymin = ymin, ymax = ymax), fill = conf_col, colour = NA) +
    geom_text(data = arcs, inherit.aes = FALSE,
              aes(x = R_arc + 0.02 * xmax, y = ymid, label = region, angle = tang, hjust = hj),
              size = 2.4, fontface = "italic", colour = "grey10")
}
p <- p +
  labs(caption = paste0("Red branches/arcs: splits absent from the genome tree. ", comparison_note)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.42, "cm"),
        plot.caption = element_text(size = 8, hjust = 0))

## ---- 10. panel b: one violin plot, three subphyla -----------------------------
sub_levels <- c("Mucoromycotina", "Mortierellomycotina", "Glomeromycotina")
dd <- meta2 %>% filter(role == "ingroup") %>%
  mutate(group = factor(group, levels = sub_levels), value = .data[[violin_var]])
nlab <- dd %>% count(group)
xlabs <- setNames(paste0(nlab$group, "\n(n = ", nlab$n, ")"), as.character(nlab$group))
pb <- ggplot(dd, aes(group, value, fill = group, colour = group)) +
  geom_violin(width = 0.8, alpha = 0.25, colour = NA, scale = "width") +
  geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.35) +
  geom_jitter(width = 0.08, height = 0, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_colour_manual(values = grp_cols, guide = "none") +
  scale_x_discrete(labels = xlabs) +
  labs(x = NULL, y = violin_ylab) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(size = 7, face = "italic"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.tag = element_text(size = 13, face = "bold"))

## ---- 11. save -----------------------------------------------------------------
tagt <- theme(plot.tag = element_text(size = 13, face = "bold"))
if (panel_b_position == "inset") {
  full <- (p + labs(tag = "a") + tagt) +
    inset_element(pb + labs(tag = "b"), left = 0.80, bottom = 0.02, right = 1.00, top = 0.30, align_to = "full")
  ggsave(paste0(out_prefix, ".pdf"), full, width = fig_width_in, height = fig_height_in, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), full, width = fig_width_in, height = fig_height_in, dpi = 400, limitsize = FALSE)
} else {
  full <- (p + labs(tag = "a") + tagt) / (plot_spacer() | (pb + labs(tag = "b")) | plot_spacer()) +
    plot_layout(heights = c(3.6, 1))
  ggsave(paste0(out_prefix, ".pdf"), full, width = fig_width_in, height = fig_height_in + 3.5, limitsize = FALSE)
  ggsave(paste0(out_prefix, ".png"), full, width = fig_width_in, height = fig_height_in + 3.5, dpi = 400, limitsize = FALSE)
}
ggsave(paste0(out_prefix, "_a_tree.pdf"), p, width = fig_width_in, height = fig_height_in, limitsize = FALSE)
ggsave(paste0(out_prefix, "_b_violin.pdf"), pb, width = 3.6, height = 3.2)
med <- dd %>% group_by(group) %>% summarise(n = n(), median = median(value), q1 = quantile(value, .25), q3 = quantile(value, .75))
print(med)
cat("\n==== 完成 ====  标记的差异分支:", length(conf_nodes), "（应为 7）\n")
