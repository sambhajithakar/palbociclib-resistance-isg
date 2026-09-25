import numpy as np, pandas as pd, networkx as nx, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig_style import *
O = "/home/claude/an/out"
E = pd.read_csv(f"{O}/ORA_results.csv"); G = pd.read_csv(f"{O}/GSEA_Hallmark.csv")
C = pd.read_csv(f"{O}/network_centrality.csv", index_col=0); rob = pd.read_csv(f"{O}/robust_DEGs.csv")

fig = plt.figure(figsize=(7.2, 8.0))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.5], hspace=0.8, wspace=0.55)

def barh(ax, d, col, title, letter, dx=-0.55, n=8):
    d = d.nsmallest(n, "p").iloc[::-1]
    lab = [t.split(" (GO:")[0][:42] for t in d.term]
    y = np.arange(len(d))
    ax.barh(y, -np.log10(d.p), color=col, height=0.62)
    ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=5.8)
    for i, (v, k, s) in enumerate(zip(-np.log10(d.p), d.overlap, d.set_size)):
        ax.text(v + 0.15, i, f"{k}/{s}", va="center", fontsize=5.2, color=INK2)
    ax.set_xlabel("$-$log$_{10}$ $P$"); ax.set_title(title, fontsize=7.2)
    ax.set_xlim(0, (-np.log10(d.p)).max() * 1.22); panel(ax, letter, dx=dx)

barh(fig.add_subplot(gs[0, 0]), E[(E.direction == "Up") & (E.library == "Hallmark")], UP,
     "Hallmark — genes up in resistance", "A")
barh(fig.add_subplot(gs[0, 1]), E[(E.direction == "Up") & (E.library == "GO_BP")], UP,
     "GO biological process — up", "B", dx=-0.62)
barh(fig.add_subplot(gs[1, 0]), E[(E.direction == "Up") & (E.library == "KEGG")], UP,
     "KEGG — up", "C")
# D GSEA NES
axd = fig.add_subplot(gs[1, 1])
g = G[G.FDR < 0.05].copy()
g = pd.concat([g.nlargest(7, "NES"), g.nsmallest(5, "NES")]).sort_values("NES")
y = np.arange(len(g))
axd.barh(y, g.NES, color=[UP if v > 0 else DOWN for v in g.NES], height=0.62)
axd.set_yticks(y); axd.set_yticklabels([t[:34] for t in g.term], fontsize=5.8)
axd.axvline(0, color=INK2, lw=0.6)
axd.set_xlabel("Normalised enrichment score"); axd.set_title("Hallmark GSEA (FDR < 0.05)", fontsize=7.2)
axd.legend(handles=[Line2D([], [], color=UP, lw=4, label="Enriched in resistant"),
                    Line2D([], [], color=DOWN, lw=4, label="Depleted in resistant")],
           fontsize=5.5, loc="lower right")
panel(axd, "D", dx=-0.62)

# E network (hub-centred subnetwork) + F hub centrality
sub = fig.add_gridspec(1, 2, width_ratios=[1.9, 1], wspace=0.42,
                       left=gs[2, :].get_position(fig).x0, right=gs[2, :].get_position(fig).x1,
                       bottom=gs[2, :].get_position(fig).y0, top=gs[2, :].get_position(fig).y1)
axe = fig.add_subplot(sub[0, 0])
g2 = nx.read_gml(f"{O}/ppi_lcc.gml")
hubs = list(C.index[C.n_top >= 3])
keep = set(hubs)
for h in hubs: keep |= set(g2.neighbors(h))
gs_ = g2.subgraph(keep).copy()
pos = nx.kamada_kawai_layout(gs_)
dirn = rob.set_index("gene").direction
deg_full = dict(g2.degree())
nx.draw_networkx_edges(gs_, pos, ax=axe, edge_color="#DCDCDC", width=0.35)
cols = [UP if dirn.get(n) == "Up" else DOWN for n in gs_.nodes()]
sizes = [14 + 3.2 * deg_full[n] for n in gs_.nodes()]
edgec = ["#111111" if n in hubs else "white" for n in gs_.nodes()]
lw = [0.9 if n in hubs else 0.25 for n in gs_.nodes()]
nx.draw_networkx_nodes(gs_, pos, ax=axe, node_color=cols, node_size=sizes, linewidths=lw, edgecolors=edgec)
lab_hubs = list(C.index[C.n_top == 4])
P = np.array([pos[n] for n in lab_hubs]); ctr = np.array([pos[n] for n in gs_.nodes()]).mean(0)
for n, p in zip(lab_hubs, P):
    v = p - ctr; v = v / (np.linalg.norm(v) + 1e-9)
    axe.annotate(n, p, fontsize=5.6, fontweight="bold", ha="center", va="center",
                 xytext=(15 * v[0], 15 * v[1]), textcoords="offset points",
                 bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8),
                 arrowprops=dict(arrowstyle="-", lw=0.4, color="#999999", shrinkA=0, shrinkB=2))
axe.set_axis_off(); axe.set_aspect("equal")
axe.set_title(f"Hub subnetwork ({gs_.number_of_nodes()} of {g2.number_of_nodes()} nodes)\n(hubs ranked top-20 by all four metrics are labelled)", fontsize=7.2)
axe.legend(handles=[Line2D([], [], marker="o", color="w", markerfacecolor=UP, ms=5, label="Up"),
                    Line2D([], [], marker="o", color="w", markerfacecolor=DOWN, ms=5, label="Down"),
                    Line2D([], [], marker="o", color="w", markerfacecolor="#CCCCCC",
                           markeredgecolor="#111111", ms=6, label="Hub")],
           fontsize=5.8, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06))
panel(axe, "E", dx=0.0, dy=1.04)

axf = fig.add_subplot(sub[0, 1])
h = C.loc[hubs].sort_values("Degree")
y = np.arange(len(h))
axf.barh(y, h.Degree, color=UP, height=0.62)
axf.set_yticks(y); axf.set_yticklabels(h.index, fontsize=5.8)
for i, v in enumerate(h.Degree): axf.text(v + 0.6, i, int(v), va="center", fontsize=5.2, color=INK2)
axf.set_xlabel("Degree in STRING network"); axf.set_title("Hub genes", fontsize=7.2)
axf.set_xlim(0, h.Degree.max() * 1.18)
panel(axf, "F", dx=-0.52, dy=1.04)
save(fig, "Fig3_function")
print(C.head(16).round(3).to_string())
