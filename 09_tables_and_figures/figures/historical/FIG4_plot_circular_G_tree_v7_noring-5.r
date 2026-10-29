## =============================================================================
## Circular species tree - G tree (genome-based), v7 style, NO BRAKER ring
## - same style as plot_circular_species_tree_v7.r (P tree)
## - tree + species names form ONE continuous coloured sector per lineage
## - names: bold; genus + epithet italic, "sp." and strain upright; aligned ends
## - BRAKER2/BRAKER3 annotation ring REMOVED (G tree is genome-based)
## 用法: RStudio 打开本脚本 -> 点 "Source"。
##       若脚本所在文件夹里找不到 G_*.treefile / tip_metadata*.tsv，会弹窗让你选。
## =============================================================================

## ---- 0. settings you may change --------------------------------------------
tree_file  <- ""          # 留空=自动找/弹窗；也可直接写路径，如 "~/Desktop/G_occ80_g90.treefile"
meta_file  <- ""          # 留空=自动找/弹窗；即 P 树用的 tip_metadata.tsv
label_size       <- 1.5
shade_alpha      <- 0.14  # lineage shading (tree + names), same value = seamless
label_band_scale <- 1.00  # >1 if names overflow, <1 if too much empty space
keep_outgroups   <- TRUE
open_angle       <- 12
fig_width_in     <- 14
fig_height_in    <- 13
weak_shalrt      <- 80    # SH-aLRT < 80  -> 低支持
weak_ufboot      <- 95    # UFBoot  < 95  -> 低支持
out_prefix       <- ""    # 留空 = "Fig_<树文件名>_v7_noring"

`%||%` <- function(a, b) if (length(a) && !is.na(a[1]) && nzchar(a[1])) a else b

## ---- 1. packages -----------------------------------------------------------
pkgs_cran <- c("ape", "ggplot2", "dplyr", "readr", "stringr", "rstudioapi")
pkgs_bioc <- c("ggtree", "treeio")
need <- setdiff(c(pkgs_cran, pkgs_bioc), rownames(installed.packages()))
if (length(need)) {
  if (length(intersect(need, pkgs_cran))) install.packages(intersect(need, pkgs_cran))
  if (length(intersect(need, pkgs_bioc))) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(intersect(need, pkgs_bioc), update = FALSE, ask = FALSE)
  }
}
suppressPackageStartupMessages({
  library(ape); library(treeio); library(ggtree)
  library(ggplot2); library(dplyr); library(readr); library(stringr)
})

## ---- 2. locate files -------------------------------------------------------
script_dir <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                       error = function(e) "")
search_dirs <- unique(path.expand(c(script_dir,
  "~/Desktop/Dual_tree_occ80_g90", "~/桌面/Dual_tree_occ80_g90", getwd())))
search_dirs <- search_dirs[nzchar(search_dirs) & dir.exists(search_dirs)]

if (!nzchar(tree_file)) {
  hits <- unlist(lapply(search_dirs, list.files, pattern = "^G_.*treefile$", full.names = TRUE))
  if (length(hits)) tree_file <- hits[1] else {
    message("请在弹出窗口里选中 G 树的 treefile 文件")
    tree_file <- file.choose()
  }
}
tree_file <- normalizePath(path.expand(tree_file))
work_dir  <- dirname(tree_file)

if (!nzchar(meta_file)) {
  hits <- unlist(lapply(unique(c(work_dir, search_dirs)), list.files,
                        pattern = "tip_metadata.*\\.tsv$", full.names = TRUE))
  if (length(hits)) meta_file <- hits[1] else {
    message("请在弹出窗口里选中 tip_metadata.tsv")
    meta_file <- file.choose()
  }
}
meta_file <- normalizePath(path.expand(meta_file))
setwd(work_dir)
if (!nzchar(out_prefix))
  out_prefix <- paste0("Fig_", sub("\\.treefile$|\\.tre$|\\.nwk$", "", basename(tree_file)), "_v7_noring")
message("工作文件夹: ", work_dir, "\n树文件: ", basename(tree_file),
        "\n信息表: ", basename(meta_file), "\n输出前缀: ", out_prefix)

## ---- 3. read, check, root --------------------------------------------------
tr <- read.tree(tree_file)
tr$tip.label <- gsub("^'|'$", "", trimws(tr$tip.label))
meta <- read_tsv(meta_file, show_col_types = FALSE, na = c("", "NA"),
                 col_types = cols(.default = col_character()))

bad1 <- setdiff(tr$tip.label, meta$label)          # 树有表无 -> 必须报错
bad2 <- setdiff(meta$label, tr$tip.label)          # 表有树无 -> 只提示
if (length(bad1)) stop("树里这些名字在信息表中找不到: ", paste(bad1, collapse = ", "))
if (length(bad2)) message("提示: 信息表中 ", length(bad2), " 个名字不在 G 树中(忽略): ",
                          paste(bad2, collapse = ", "))
meta <- meta %>% filter(label %in% tr$tip.label)

og <- meta$label[meta$role == "outgroup"]
if (length(og) && is.monophyletic(tr, og)) {
  tr <- root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
} else if (length(og) && is.monophyletic(tr, setdiff(tr$tip.label, og))) {
  tr <- root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
} else if (length(og)) {
  asc <- intersect(meta$label[meta$outer_phylum %in% c("Ascomycota", "Basidiomycota")], og)
  warning("外部参照在 G 树中不成一支，改用 Dikarya 外类群定根: ", paste(asc, collapse = ", "))
  tr <- root(tr, outgroup = asc, resolve.root = TRUE, edgelabel = TRUE)
}
if (!keep_outgroups && length(og)) tr <- drop.tip(tr, og)
tr <- ladderize(tr, right = FALSE)

## 支持度: 兼容 "SH-aLRT/UFBoot" 或 单个数值(视为 UFBoot)
nl <- if (is.null(tr$node.label)) rep("", tr$Nnode) else tr$node.label
sp <- str_split_fixed(nl, "/", 2)
v1 <- suppressWarnings(as.numeric(sp[, 1])); v2 <- suppressWarnings(as.numeric(sp[, 2]))
two <- !is.na(v2)
sup <- tibble(node = Ntip(tr) + seq_along(nl),
              shalrt = ifelse(two, v1, NA), ufboot = ifelse(two, v2, v1)) %>%
  mutate(weak = (!is.na(shalrt) & shalrt < weak_shalrt) | (!is.na(ufboot) & ufboot < weak_ufboot))
weak_nodes <- sup$node[sup$weak]

tr_plot <- compute.brlen(tr, method = "Grafen", power = 0.7)
tr_plot$node.label <- tr$node.label

## ---- 4. scientific-name formatting ----------------------------------------
fmt_name <- function(x) {
  x <- gsub("\\s+", " ", trimws(x))
  x <- gsub("\\b(sp|cf|aff)\\.(?=\\S)", "\\1. ", x, perl = TRUE)
  esc <- function(s) gsub('"', '\\\\"', s)
  m <- regmatches(x, regexec("^(\\S+)\\s+(sp\\.|cf\\.|aff\\.)\\s*(.*)$", x))[[1]]
  if (length(m)) {
    sprintf('bolditalic("%s")~bold("%s")', esc(m[2]), esc(trimws(paste(m[3], m[4]))))
  } else sprintf('bolditalic("%s")', esc(x))
}

meta2 <- meta %>%
  filter(label %in% tr_plot$tip.label) %>%
  mutate(display_name = ifelse(is.na(display_name), gsub("_", " ", label), display_name),
         group = ifelse(is.na(group) & role == "outgroup", "Outgroup", group),
         name_expr = vapply(display_name, fmt_name, character(1)))

## ---- 5. colours ------------------------------------------------------------
grp_levels <- c("Glomeromycotina", "Mortierellomycotina", "Mucoromycotina", "Outgroup")
grp_cols <- c(Glomeromycotina = "#4DAF4A", Mortierellomycotina = "#984EA3",
              Mucoromycotina = "#377EB8", Outgroup = "grey55")
grp_cols <- grp_cols[names(grp_cols) %in% unique(meta2$group)]

hl <- bind_rows(lapply(intersect(grp_levels, unique(meta2$group)), function(g) {
  tp <- meta2$label[meta2$group == g]
  if (length(tp) == 1) return(data.frame(lineage = g, node = which(tr_plot$tip.label == tp)))
  if (!is.monophyletic(tr_plot, tp, reroot = FALSE)) {
    message("注意: ", g, " 在 G 树中不是单系，树内部不涂色，只给名字带涂色")
    return(NULL)
  }
  data.frame(lineage = g, node = getMRCA(tr_plot, tp))
}))

## ---- 6. geometry -----------------------------------------------------------
p0 <- ggtree(tr_plot, layout = "fan", open.angle = open_angle, size = 0.3)
xmax <- max(p0$data$x[p0$data$isTip], na.rm = TRUE)

has_order <- "order" %in% names(meta2) && any(!is.na(meta2$order))
outer_pad <- 0.02 + (if (has_order) 0.12 else 0)

max_chars  <- max(nchar(meta2$display_name))
label_in   <- max_chars * label_size * 0.62 / 25.4 * label_band_scale
panel_r_in <- 0.46 * min(fig_width_in - 2.8, fig_height_in)
ratio  <- min(0.85, label_in * (1 + outer_pad) / panel_r_in)
R_in   <- xmax                      # names start exactly where the tree ends
R_out  <- R_in / (1 - ratio)
R_text <- R_out - 0.008 * R_out

tips <- p0$data %>% filter(isTip) %>% select(label, y, angle) %>%
  left_join(meta2 %>% select(label, name_expr, group), by = "label") %>%
  mutate(flip = angle > 90 & angle < 270,
         tangle = ifelse(flip, angle + 180, angle),
         hj = ifelse(flip, 0, 1))

## ---- 7. tree + seamless lineage sector (tree shading + name band) ----------
p <- p0
if (nrow(hl)) p <- p + geom_hilight(data = hl, mapping = aes(node = node, fill = lineage),
                                    alpha = shade_alpha)
p <- p +
  geom_rect(data = tips, inherit.aes = FALSE,
            aes(xmin = R_in, xmax = R_out, ymin = y - 0.5, ymax = y + 0.5, fill = group),
            alpha = shade_alpha, colour = NA) +
  scale_fill_manual(values = grp_cols, breaks = names(grp_cols), name = "Lineage") +
  geom_point2(aes(subset = node %in% weak_nodes),
              shape = 21, size = 1.4, fill = "white", colour = "black", stroke = 0.4) +
  geom_text(data = tips, inherit.aes = FALSE, parse = TRUE,
            aes(x = R_text, y = y, label = name_expr, angle = tangle, hjust = hj),
            size = label_size, colour = "grey10")

## ---- 8. optional order arcs (only if metadata has an "order" column) -------
if (has_order) {
  ord <- split(meta2$label[!is.na(meta2$order)], meta2$order[!is.na(meta2$order)])
  clade_df <- bind_rows(lapply(names(ord), function(o) {
    tp <- ord[[o]]
    node <- if (length(tp) == 1) which(tr_plot$tip.label == tp) else getMRCA(tr_plot, tp)
    mono <- length(tp) == 1 || is.monophyletic(tr_plot, tp, reroot = FALSE)
    data.frame(node = node, clade_name = o, n = length(tp), mono = mono)
  }))
  if (any(!clade_df$mono))
    message("以下目在树中不是单系，未画弧线: ", paste(clade_df$clade_name[!clade_df$mono], collapse = ", "))
  clade_df <- clade_df[clade_df$mono, ]
  p <- p + geom_cladelab(data = clade_df, mapping = aes(node = node, label = clade_name),
                         offset = R_out - xmax + 0.02 * xmax, offset.text = 0.02 * xmax,
                         barsize = 0.8, fontsize = 2.6, angle = "auto", hjust = 0.5,
                         horizontal = FALSE, barcolour = "grey30", textcolour = "grey15")
}

p <- p + theme(legend.position = "right",
               legend.title = element_text(size = 9, face = "bold"),
               legend.text  = element_text(size = 8),
               legend.key.size = unit(0.42, "cm"))

## ---- 9. save ---------------------------------------------------------------
ggsave(paste0(out_prefix, ".pdf"), p, width = fig_width_in, height = fig_height_in, limitsize = FALSE)
ggsave(paste0(out_prefix, ".png"), p, width = fig_width_in, height = fig_height_in, dpi = 400, limitsize = FALSE)
print(p)
cat("\n==== 完成 ====\n已保存到", work_dir, ":", paste0(out_prefix, c(".pdf", ".png")), "\n")
cat("叶片数:", Ntip(tr_plot), " | 低支持分支(空心圆):", length(weak_nodes), "\n")
