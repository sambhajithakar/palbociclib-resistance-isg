import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from fig_style import *

fig, ax = plt.subplots(figsize=(7.2, 4.6)); ax.set_xlim(0, 100); ax.set_ylim(0, 72); ax.axis("off")
FILL = {"data": "#EAF2F8", "an": "#FDF0E6", "val": "#E8F5F0", "out": "#F2EDF6"}
EDGE = {"data": CAT[0], "an": CAT[1], "val": CAT[2], "out": CAT[3]}

def box(x, y, w, h, title, body, kind):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.4",
                                fc=FILL[kind], ec=EDGE[kind], lw=1.0))
    ax.text(x + w / 2, y + h - 2.6, title, ha="center", va="top", fontsize=7.2, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - 7.2, body, ha="center", va="top", fontsize=6.0, color=INK2, linespacing=1.5)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=9,
                                 lw=0.9, color="#8A8A8A", shrinkA=1, shrinkB=1))

box(1, 47, 30, 22, "Discovery data (GEO)",
    "GSE130437  MCF7 palbo-R\nGSE222367  MCF7 / T47D palbo-R\nGSE229235  ER+ PDX palbo-R\n49 resistant vs 25 sensitive", "data")
box(35, 47, 30, 22, "Differential expression",
    "Moderated $t$ per dataset\nRobust rank aggregation\n378 robust genes\n(264 up, 114 down)", "an")
box(69, 47, 30, 22, "Function and network",
    "GO / KEGG / Hallmark\nPreranked GSEA\nSTRING v12 network\n14 interferon hub genes", "an")

box(1, 24, 30, 19, "Machine learning",
    "LASSO · SVM-RFE · Boruta\nGradient boosting\n17-gene core signature\nNested leave-one-dataset-out", "an")
box(35, 24, 30, 19, "Clinical validation",
    "TCGA-BRCA ER+/HER2$-$ (n = 643)\nMETABRIC ER+/HER2$-$ (n = 1396)\nLASSO-Cox risk score\nNomogram · calibration · DCA", "val")
box(69, 24, 30, 19, "On-treatment biopsies",
    "NeoPalAna trial (GSE93204)\nC1D1 vs C1D15, paired\nE2F, interferon and\nresistance scores", "val")

box(18, 2, 64, 17, "Structure-based drug repurposing against STAT1",
    "STAT1 SH2 phosphotyrosine subsite (PDB 1BF5, biological dimer) · redocking validation\n"
    "1,176 approved drugs (ChEMBL) · AutoDock Vina two-stage screen · triplicate re-docking\n"
    "Drug-likeness and interaction analysis of the top-ranked compounds", "out")

for x in (31, 65): arrow(x, 58, x + 4, 58)
ax.plot([84, 84, 16], [47, 45.2, 45.2], color="#8A8A8A", lw=0.9, solid_capstyle="round")
arrow(16, 45.2, 16, 43.2)
arrow(31, 33.5, 35, 33.5); arrow(65, 33.5, 69, 33.5)
arrow(16, 24, 16, 19.4); arrow(50, 24, 50, 19.4); arrow(84, 24, 84, 19.4)
save(fig, "Fig1_workflow")
