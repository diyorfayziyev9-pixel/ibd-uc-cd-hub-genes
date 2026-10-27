# 53_fig1_workflow.py: Figure 1 (study workflow). Run: python 53_fig1_workflow.py (needs matplotlib); writes Figure1_workflow.png
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

UC, CD, CTRL, HUB, GREY = "#D1495B", "#8E3B6E", "#5B7FA6", "#C0392B", "#EEF1F5"
plt.rcParams["font.family"] = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(11, 12.5), dpi=300)
ax.set_xlim(0, 11); ax.set_ylim(0, 12.5); ax.axis("off")

def box(x, y, w, h, title, body="", edge="#444", fill="white", tcol="black", lw=1.6, tsize=11, bsize=9.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                fc=fill, ec=edge, lw=lw))
    if body:
        ax.text(x + w/2, y + h - 0.28, title, ha="center", va="top", fontsize=tsize, fontweight="bold", color=tcol)
        ax.text(x + w/2, y + h - 0.68, body, ha="center", va="top", fontsize=bsize, color="#222", linespacing=1.35)
    else:
        ax.text(x + w/2, y + h/2, title, ha="center", va="center", fontsize=tsize, fontweight="bold", color=tcol)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16, lw=1.6, color="#555"))

# 1 data
box(1.2, 10.6, 8.6, 1.5, "Discovery cohort GSE75214", "Affymetrix Human Gene 1.0 ST, biopsies\nUC n = 74 vs control n = 11     CD n = 51 vs control n = 11",
    edge="#444", fill=GREY)

# 2 DEG
box(0.6, 8.55, 4.5, 1.45, "UC vs control", "limma, adj. P < 0.05, |log2FC| > 1\n1,423 probe sets = 1,281 genes", edge=UC, tcol=UC)
box(5.9, 8.55, 4.5, 1.45, "CD vs control", "limma, adj. P < 0.05, |log2FC| > 1\n551 probe sets = 498 genes", edge=CD, tcol=CD)
arrow(3.0, 10.6, 2.85, 10.0); arrow(8.0, 10.6, 8.15, 10.0)

# 3 overlap
box(2.6, 6.6, 5.8, 1.5, "Shared DEGs: 328 genes", "321 concordant (227 up, 94 down)\n7 opposite direction (excluded)", edge="#444", fill=GREY)
arrow(2.85, 8.55, 4.3, 8.1); arrow(8.15, 8.55, 6.7, 8.1)

# 4 enrichment + PPI
box(0.6, 4.5, 4.5, 1.55, "Functional enrichment", "GO Biological Process, KEGG\n(also UC- and CD-specific sets)", edge="#444")
box(5.9, 4.5, 4.5, 1.55, "PPI network (STRING, score 0.7)", "179 connected nodes, 488 edges;\nlargest component 135 nodes, 450 edges", edge="#444")
arrow(4.3, 6.6, 2.85, 6.05); arrow(6.7, 6.6, 8.15, 6.05)

# hub
box(5.9, 2.75, 4.5, 1.2, "5 hub genes", "IL6, IL1B, CCL2, STAT1, MMP2", edge=HUB, tcol=HUB, fill="#FBEDEB")
arrow(8.15, 4.5, 8.15, 3.95)

# 5 validation / immune / drugs
box(0.2, 0.35, 3.3, 1.9, "Validation", "GSE206171 (UC)\nGSE66407 (CD ileum)\nWilcoxon + BH, ROC", edge=CTRL, tcol=CTRL)
box(3.85, 0.35, 3.3, 1.9, "Immune context", "xCell enrichment scores\nSpearman correlation\n(all, UC only, CD only)", edge="#444")
box(7.5, 0.35, 3.3, 1.9, "Drug-gene links", "DGIdb v5.0.13\nrepurposing hypotheses\n(not tested)", edge="#444")
arrow(7.2, 2.75, 1.9, 2.25); arrow(8.15, 2.75, 5.5, 2.25); arrow(9.0, 2.75, 9.15, 2.25)

fig.savefig("Figure1_workflow.png", bbox_inches="tight", facecolor="white")
print("ok")