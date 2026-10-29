# -*- coding: utf-8 -*-
"""
Fig_contamination_filtering_omark.py  (v4)
  a  OMArk 污染比例: 过滤前 (空心) -> 过滤后 (实心)
  b  OMArk 原生样式: 蛋白数 / 完整度 (S,D,M) / 一致性 (A,C,I,U; Partial/Fragment 斜线)
配色/斜线/堆叠顺序与 OMArk src/omark/graphics.py 一致
v4: 物种名不用 mathtext (PDF 中为整段可编辑文字); a-b 间距按实际字体自动计算
用法: RStudio -> reticulate::source_python("脚本路径")  或终端 python3 本文件
"""
import os, math
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib import transforms as mtrans

# ===================== 可调参数 =====================
INFILE = os.path.expanduser("~/Desktop/contamination_filtering_before_after_10strains_v2.tsv")
FIG_W_MM, FIG_H_MM = 180, 230
GAP_AB_MM = 9            # a 的横轴标题 与 b 顶部色条 之间的距离 (mm)
NAME_SIZE_A = 7          # a 中物种名字号
NAME_SIZE_B = 6.5        # b 中物种名字号
NAME_ROT_B = 60          # b 中物种名旋转角度
# ===================================================

if not os.path.exists(INFILE):
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk(); root.withdraw()
        INFILE = filedialog.askopenfilename(title="选择 contamination_filtering_before_after_10strains_v2.tsv")
        root.destroy()
    except Exception:
        pass
if not INFILE or not os.path.exists(INFILE):
    raise SystemExit("找不到输入表格, 请修改脚本顶部 INFILE")
OUTDIR = os.path.dirname(os.path.abspath(INFILE))

F_SINGLE, F_DUP, MISSING = "#46bea1", "#cfee8e", "#ed1c5a"
CONSISTENT, INCONSISTENT = "#3db7e9", "#8a5df1"
CONTAMINATION, UNKNOWN = "#e69f00", "#36454f"
PROT_B, PROT_A = "#9ec2e2", "#2f6fa7"
SUB_COLS = {"Mucoromycotina": "#1b9e77", "Mortierellomycotina": "#d95f02", "Glomeromycotina": "#7570b3"}

matplotlib.rcParams.update({"font.family": "sans-serif", "font.size": 7, "axes.linewidth": 0.6,
                            "hatch.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42,
                            "svg.fonttype": "none"})

d = pd.read_csv(INFILE, sep="\t")
need = ["sample", "species_label", "subphylum", "omark_before_C_pct", "omark_after_C_pct",
        "omark_before_A_full_pct", "omark_after_A_partial_pct", "omark_before_S_pct", "omark_before_proteins"]
miss = [c for c in need if c not in d.columns]
if miss:
    raise SystemExit(f"表格缺少列 {miss}; 请使用最新的 contamination_filtering_before_after_10strains_v2.tsv")
d["subphylum"] = pd.Categorical(d["subphylum"], list(SUB_COLS), ordered=True)
d = d.sort_values(["subphylum", "omark_before_C_pct"], ascending=[True, False]).reset_index(drop=True)

if "reclassified_flag" in d.columns:
    flag = d["reclassified_flag"].astype(str).str.lower().eq("true")
elif "note" in d.columns:
    flag = d["note"].fillna("").astype(str).str.contains("no longer detected", case=False)
else:
    flag = d["sample"].isin(["Acaulospora_colombiana", "Oehlia_diaphana"])

# ---- 物种名: (斜体段, 正体段); 每段是一个完整文字对象 ----
def split_name(s):
    if " sp. " in s:
        g, rest = s.split(" sp. ", 1)
        return g, "sp. " + rest
    return s, ""
names = [split_name(s) for s in d["species_label"]]

def draw_name(ax, x, y, italic, upright, transform, rotation=0, fontsize=7, va="center"):
    """右对齐书写 '斜体 正体'; 正体段在锚点处, 斜体段沿文字方向向左平移"""
    kw = dict(rotation=rotation, rotation_mode="anchor", fontsize=fontsize, ha="right", va=va, clip_on=False)
    if not upright:
        ax.text(x, y, italic, style="italic", transform=transform, **kw)
        return
    ax.text(x, y, upright, transform=transform, **kw)
    rnd = ax.figure.canvas.get_renderer()
    tmp = ax.figure.text(0, 0, upright, fontsize=fontsize)
    w_pt = tmp.get_window_extent(rnd).width * 72 / ax.figure.dpi + fontsize * 0.30
    tmp.remove()
    a = math.radians(rotation)
    off = mtrans.offset_copy(transform, fig=ax.figure, x=-w_pt * math.cos(a), y=-w_pt * math.sin(a), units="points")
    ax.text(x, y, italic, style="italic", transform=off, **kw)

fig = plt.figure(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4))
gsA = fig.add_gridspec(1, 1, left=0.30, right=0.74, top=0.975, bottom=0.735)
gsB = fig.add_gridspec(3, 1, height_ratios=[0.42, 1.0, 1.0], hspace=0.06,
                       left=0.11, right=0.80, top=0.62, bottom=0.165)
ax_a = fig.add_subplot(gsA[0])
ax_p = fig.add_subplot(gsB[0]); ax_c = fig.add_subplot(gsB[1], sharex=ax_p); ax_k = fig.add_subplot(gsB[2], sharex=ax_p)
fig.text(0.015, 0.985, "a", fontweight="bold", fontsize=11, va="top")
lab_b = fig.text(0.015, 0.65, "b", fontweight="bold", fontsize=11, va="top")
fig.canvas.draw()

# ===== a =====
y = list(range(len(d)))[::-1]
for yi, (_, r), f in zip(y, d.iterrows(), flag):
    ax_a.plot([r.omark_before_C_pct, r.omark_after_C_pct], [yi, yi], color="grey", lw=1, zorder=1)
    ax_a.scatter(r.omark_before_C_pct, yi, s=26, facecolor="white", edgecolor="#333333", lw=0.8, zorder=2)
    ax_a.scatter(r.omark_after_C_pct, yi, s=30, facecolor=SUB_COLS[r.subphylum], edgecolor="black", lw=0.6, zorder=3)
    if f:
        ax_a.text(r.omark_after_C_pct - 0.6, yi + 0.18, "*", ha="right", va="center", fontsize=9)
ax_a.set_yticks(y); ax_a.set_yticklabels([""] * len(y))
ax_a.set_xlim(left=-0.8); ax_a.set_ylim(-0.6, len(d) - 0.4)
ax_a.set_xlabel("OMArk contamination (% of proteins)")
ax_a.grid(axis="x", lw=0.3, color="#dddddd"); ax_a.set_axisbelow(True)
ytr = mtrans.blended_transform_factory(ax_a.transAxes, ax_a.transData)
for yi, (it, up) in zip(y, names):
    draw_name(ax_a, -0.025, yi, it, up, ytr, rotation=0, fontsize=NAME_SIZE_A, va="center")
h = [Line2D([], [], marker="o", ls="", mfc="white", mec="#333333", label="Before filtering")] + \
    [Line2D([], [], marker="o", ls="", mfc=c, mec="black", label=f"After ({k})") for k, c in SUB_COLS.items()]
ax_a.legend(handles=h, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6.5)

# ===== b =====
W, GAP = 0.9, 0.7
BW = W * 0.97
xs_b, xs_a, centers, x = [], [], [], 0.0
for _ in range(len(d)):
    xs_b.append(x); xs_a.append(x + W); centers.append(x + W); x += 2 * W + GAP
stages = [("before", xs_b), ("after", xs_a)]
FULL = dict(width=BW, align="edge", edgecolor="none", linewidth=0, alpha=.99)

def hbar(ax, xs, h, bottom, color, kind):
    if kind == "full":
        ax.bar(xs, h, bottom=bottom, color=color, **FULL)
    else:
        ax.bar(xs, h, bottom=bottom, facecolor=color, hatch="///",
               edgecolor=("darkslategray" if kind == "partial" else "whitesmoke"),
               linewidth=0.1, width=BW, align="edge", alpha=.99)
    return bottom + h

for st, xs in stages:
    ax_p.bar(xs, d[f"omark_{st}_proteins"] / 1000, color=PROT_A if st == "after" else PROT_B, **FULL)
ax_p.set_ylabel("Proteins\n(×1,000)")

for st, xs in stages:
    S, D, M = d[f"omark_{st}_S_pct"], d[f"omark_{st}_D_pct"], d[f"omark_{st}_M_pct"]
    ax_c.bar(xs, S, color=F_SINGLE, **FULL)
    ax_c.bar(xs, D, bottom=S, color=F_DUP, **FULL)
    ax_c.bar(xs, M, bottom=S + D, color=MISSING, **FULL)
ax_c.set_ylim(0, 100); ax_c.set_ylabel("Percentage of\nconserved HOGs")

for st, xs in stages:
    btm = pd.Series(0.0, index=d.index)
    for cat, col in (("A", CONSISTENT), ("C", CONTAMINATION), ("I", INCONSISTENT)):
        for kind in ("full", "partial", "fragment"):
            btm = hbar(ax_k, xs, d[f"omark_{st}_{cat}_{kind}_pct"], btm, col, kind)
    ax_k.bar(xs, (100 - btm).clip(lower=0), bottom=btm, color=UNKNOWN, **FULL)
ax_k.set_ylim(100, 0); ax_k.set_ylabel("Percentage of\nproteome")

for ax in (ax_p, ax_c):
    plt.setp(ax.get_xticklabels(), visible=False); ax.tick_params(axis="x", length=0)
for ax in (ax_p, ax_c, ax_k):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
ax_c.spines["bottom"].set_linewidth(1.0)
ax_k.set_xticks([xx + BW / 2 for xx in xs_b + xs_a])
ax_k.set_xticklabels(["B"] * len(d) + [("A*" if f else "A") for f in flag], fontsize=5.5)
ax_k.tick_params(axis="x", length=0, pad=1.5)
ax_k.set_xlim(-GAP / 2, x - GAP / 2)

STRIP = 0.12   # 亚门色条高度 (占蛋白面板高度比例)
for c, sub in zip(centers, d["subphylum"]):
    ax_p.add_patch(mpatches.Rectangle((c - W, 1.05), 2 * W - 0.02, STRIP * 0.45, transform=ax_p.get_xaxis_transform(),
                                      color=SUB_COLS[sub], clip_on=False))

ax_p.legend(handles=[mpatches.Patch(color=PROT_B, label="Before (B)"), mpatches.Patch(color=PROT_A, label="After (A)")],
            title="Proteins", fontsize=6, title_fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1.0))
ax_c.legend(handles=[mpatches.Patch(color=F_SINGLE, label="Single"), mpatches.Patch(color=F_DUP, label="Duplicated"),
                     mpatches.Patch(color=MISSING, label="Missing")],
            title="Completeness", fontsize=6, title_fontsize=6.5, loc="center left", bbox_to_anchor=(1.01, 0.5))
ax_k.legend(handles=[mpatches.Patch(color=CONSISTENT, label="Consistent"),
                     mpatches.Patch(color=CONTAMINATION, label="Contamination"),
                     mpatches.Patch(color=INCONSISTENT, label="Inconsistent"),
                     mpatches.Patch(color=UNKNOWN, label="Unknown"),
                     mpatches.Patch(facecolor="whitesmoke", hatch="///", edgecolor="darkslategray", lw=0.1, label="Partial mapping"),
                     mpatches.Patch(facecolor="darkslategray", hatch="///", edgecolor="whitesmoke", lw=0.1, label="Fragments")],
            title="Consistency", fontsize=6, title_fontsize=6.5, loc="center left", bbox_to_anchor=(1.01, 0.5))

# ---- 按实际字体测量 a 的下边界, 确定 b 顶边 ----
fig.canvas.draw()
rnd = fig.canvas.get_renderer()
a_bottom = ax_a.get_tightbbox(rnd).transformed(fig.transFigure.inverted()).y0
H = FIG_H_MM
p_h = ax_p.get_position().height
new_top = a_bottom - GAP_AB_MM / H - p_h * 0.12
gsB.update(top=new_top)
fig.canvas.draw()
p_h = ax_p.get_position().height
new_top = a_bottom - GAP_AB_MM / H - p_h * 0.12
gsB.update(top=new_top)
lab_b.set_y(new_top + p_h * 0.12 + 0.55 * GAP_AB_MM / H)

# ---- b 物种名 (在最终位置上绘制) ----
for c, (it, up) in zip(centers, names):
    draw_name(ax_k, c, 105, it, up, ax_k.transData, rotation=NAME_ROT_B, fontsize=NAME_SIZE_B, va="top")

fig.text(0.015, 0.012, "* Contaminant lineage no longer detected after filtering; residual proteins were reassigned "
         "to other categories by OMArk (see Supplementary Table).", fontsize=5.5, style="italic")

out_pdf = os.path.join(OUTDIR, "Fig_contamination_filtering_omark.pdf")
out_png = os.path.join(OUTDIR, "Fig_contamination_filtering_omark.png")
out_svg = os.path.join(OUTDIR, "Fig_contamination_filtering_omark.svg")
fig.savefig(out_pdf); fig.savefig(out_svg); fig.savefig(out_png, dpi=300)
print("读取表格:", INFILE)
for f in (out_pdf, out_svg, out_png): print("已保存:", f)
plt.show()
