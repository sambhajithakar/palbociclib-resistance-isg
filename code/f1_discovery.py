import numpy as np, pandas as pd, matplotlib.pyplot as plt, json
from matplotlib.lines import Line2D
from fig_style import *
O = "/home/claude/an/out"
DS = ["GSE130437", "GSE222367", "GSE229235"]
LAB = {"GSE130437": "GSE130437 (MCF7)", "GSE222367": "GSE222367 (MCF7/T47D)", "GSE229235": "GSE229235 (PDX)"}

deg = {d: pd.read_csv(f"{O}/DEG_{d}.csv") for d in DS}
rob = pd.read_csv(f"{O}/robust_DEGs.csv"); ra = pd.read_csv(f"{O}/RRA_all.csv")
M = pd.read_pickle(f"{O}/merged_centered.pkl"); meta = pd.read_csv(f"{O}/meta.csv").set_index("sample").loc[M.columns]

fig = plt.figure(figsize=(7.2, 7.6))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1.15, 1], hspace=0.62, wspace=0.58)

# A-C volcanoes
for i, d in enumerate(DS):
    ax = fig.add_subplot(gs[0, i]); r = deg[d].dropna(subset=["FDR"])
    ns = r[r.sig == "NS"]; up = r[r.sig == "Up"]; dn = r[r.sig == "Down"]
    ax.scatter(ns.log2FC, -np.log10(ns.P), s=1.2, c=NEUTRAL, lw=0, rasterized=True)
    ax.scatter(dn.log2FC, -np.log10(dn.P), s=1.6, c=DOWN, lw=0, rasterized=True)
    ax.scatter(up.log2FC, -np.log10(up.P), s=1.6, c=UP, lw=0, rasterized=True)
    for x in (-1, 1): ax.axvline(x, ls=(0, (3, 3)), c=MUTED, lw=0.5)
    ax.set_title(LAB[d], fontsize=7); ax.set_xlabel("log$_2$ fold change"); ax.set_xlim(-9, 9)
    if i == 0: ax.set_ylabel("$-$log$_{10}$ $P$")
    ax.text(0.03, 0.97, f"{len(up)} up", transform=ax.transAxes, color=UP, va="top", fontsize=6)
    ax.text(0.03, 0.88, f"{len(dn)} down", transform=ax.transAxes, color=DOWN, va="top", fontsize=6)
    panel(ax, "ABC"[i], dx=-0.28 if i == 0 else -0.2)

# D heatmap of top RRA genes
axd = fig.add_subplot(gs[1, :2])
top = pd.concat([rob[rob.direction == "Up"].nsmallest(16, "score"), rob[rob.direction == "Down"].nsmallest(10, "score")])
mat = top.set_index("gene")[[f"LFC_{d}" for d in DS]]
im = axd.imshow(mat.T.values, cmap="RdBu_r", vmin=-5, vmax=5, aspect="auto")
axd.set_xticks(range(len(mat))); axd.set_xticklabels(mat.index, rotation=90, fontsize=5.4)
axd.set_yticks(range(3)); axd.set_yticklabels([LAB[d].split(" (")[0] for d in DS], fontsize=6)
axd.set_title("Top robust genes (rank aggregation)", fontsize=7.5)
for s in axd.spines.values(): s.set_visible(True)
cb = fig.colorbar(im, ax=axd, fraction=0.025, pad=0.04); cb.set_label("log$_2$ FC", fontsize=6, labelpad=2); cb.ax.tick_params(labelsize=5.5)
panel(axd, "D", dx=-0.09, dy=1.22)

# E direction summary
axe = fig.add_subplot(gs[2, 2])
cnt = [(rob.direction == "Up").sum(), (rob.direction == "Down").sum()]
b = axe.bar(["Up", "Down"], cnt, color=[UP, DOWN], width=0.55)
for r, v in zip(b, cnt): axe.text(r.get_x() + r.get_width() / 2, v + 4, str(v), ha="center", fontsize=6.5)
axe.set_ylabel("Robust genes", labelpad=1); axe.set_title(f"{len(rob)} robust genes", fontsize=7.5); axe.set_ylim(0, max(cnt) * 1.2)
panel(axe, "H", dx=-0.42)

# F,G PCA before and after within-model centring
expr = pd.read_pickle(f"{O}/expr_all.pkl")
uni = M.index
raw = pd.concat([expr[d].reindex(uni)[meta.index[meta.dataset == d]] for d in DS], axis=1)
mk = {"Sensitive": "o", "Resistant": "^"}
def pca_panel(ax, mat, title, letter, dx):
    v = mat.loc[mat.var(1).nlargest(2000).index]
    U, S, _ = np.linalg.svd((v.T - v.mean(1)).values, full_matrices=False)
    sc = U[:, :2] * S[:2]; ev = S ** 2 / (S ** 2).sum() * 100
    mm = meta.loc[mat.columns]
    for j, d in enumerate(DS):
        for g in mk:
            k = (mm.dataset == d).values & (mm.group == g).values
            ax.scatter(sc[k, 0], sc[k, 1], s=20, c=CAT[j], marker=mk[g], lw=0.4, edgecolor="white")
    ax.set_xlabel(f"PC1 ({ev[0]:.0f}%)"); ax.set_ylabel(f"PC2 ({ev[1]:.0f}%)"); ax.set_title(title, fontsize=7.5)
    panel(ax, letter, dx=dx)
axf = fig.add_subplot(gs[2, 0]); pca_panel(axf, raw, "Before centring", "F", -0.30)
axg2 = fig.add_subplot(gs[2, 1]); pca_panel(axg2, M, "After within-model centring", "G", -0.22)
h = [Line2D([], [], color=CAT[j], marker="s", ls="", ms=5, label=d.replace("GSE", "GSE ")) for j, d in enumerate(DS)]
h += [Line2D([], [], color=MUTED, marker=mk[g], ls="", ms=5, label=g) for g in mk]
axg2.legend(handles=h, ncol=2, loc="upper center", fontsize=5.4, handletextpad=0.25, columnspacing=0.7, labelspacing=0.25)

# G sample counts
axg = fig.add_subplot(gs[1, 2])
t = meta.groupby(["dataset", "group"]).size().unstack(fill_value=0)[["Sensitive", "Resistant"]]
y = np.arange(len(t))
axg.barh(y - 0.17, t.Sensitive, 0.32, color=LOW, label="Sensitive")
axg.barh(y + 0.17, t.Resistant, 0.32, color=HIGH, label="Resistant")
for i in y:
    axg.text(t.Sensitive.iloc[i] + 0.4, i - 0.17, int(t.Sensitive.iloc[i]), va="center", fontsize=6)
    axg.text(t.Resistant.iloc[i] + 0.4, i + 0.17, int(t.Resistant.iloc[i]), va="center", fontsize=6)
axg.set_yticks(y); axg.set_yticklabels([s.replace("GSE", "") for s in t.index], fontsize=6)
axg.set_xlabel("Samples"); axg.set_title("Discovery samples", fontsize=7.5); axg.legend(fontsize=6, loc="lower right")
axg.set_xlim(0, t.values.max() * 1.25)
panel(axg, "E", dx=-0.42, dy=1.22)
save(fig, "Fig2_discovery")
print(mat.round(2).to_string())
