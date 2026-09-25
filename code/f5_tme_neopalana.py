import numpy as np, pandas as pd, json, matplotlib.pyplot as plt
from scipy import stats
from matplotlib.lines import Line2D
from fig_style import *
O = "/home/claude/an/out"
hw = pd.read_csv(f"{O}/tcga_hallmark_vs_risk.csv"); ck = pd.read_csv(f"{O}/tcga_checkpoints.csv")
T = pd.read_csv(f"{O}/tcga_cohort.csv", index_col=0)
np_sc = pd.read_csv(f"{O}/neopalana_scores.csv"); npr = json.load(open(f"{O}/neopalana_results.json"))

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.55, height_ratios=[1.25, 1])

# A Hallmark vs risk
axa = fig.add_subplot(gs[0, 0])
s = hw[hw.FDR < 0.05].copy()
s = pd.concat([s.nsmallest(8, "rho"), s.nlargest(8, "rho")]).sort_values("rho")
y = np.arange(len(s))
axa.barh(y, s.rho, color=[UP if v > 0 else DOWN for v in s.rho], height=0.66)
axa.set_yticks(y); axa.set_yticklabels([t[:30] for t in s.pathway], fontsize=5.6)
axa.axvline(0, color=INK2, lw=0.6)
axa.set_xlabel("Spearman $\\rho$ with risk score"); axa.set_title("Hallmark activity vs risk score\n(TCGA ER+/HER2$-$)", fontsize=7.2)
panel(axa, "A", dx=-0.72)

# B checkpoints
axb = fig.add_subplot(gs[0, 1:])
ck = ck.sort_values("median_diff_high_minus_low")
y = np.arange(len(ck))
cols = [UP if (d > 0 and p < 0.05) else DOWN if (d < 0 and p < 0.05) else "#C4C4C4"
        for d, p in zip(ck.median_diff_high_minus_low, ck.FDR)]
axb.barh(y, ck.median_diff_high_minus_low, color=cols, height=0.64)
axb.set_yticks(y); axb.set_yticklabels(ck.gene, fontsize=6)
axb.axvline(0, color=INK2, lw=0.6)
for i, (v, p) in enumerate(zip(ck.median_diff_high_minus_low, ck.FDR)):
    st = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    axb.text(v + (0.03 if v > 0 else -0.03), i, st, va="center", ha="left" if v > 0 else "right", fontsize=5.6)
axb.set_xlabel("Median difference in log$_2$ expression (high $-$ low risk)")
axb.set_title("Immune-checkpoint and cytotoxicity genes by risk group", fontsize=7.2)
axb.set_xlim(-1.25, 0.75)
axb.legend(handles=[Line2D([], [], color=UP, lw=4, label="Higher in high risk"),
                    Line2D([], [], color=DOWN, lw=4, label="Lower in high risk"),
                    Line2D([], [], color="#C4C4C4", lw=4, label="Not significant")],
           fontsize=5.8, loc="lower right")
panel(axb, "B", dx=-0.2)

# C-E NeoPalAna paired changes
w = np_sc.pivot_table(index="patient", columns="tp", values=["e2f_score", "resist_score", "ifn_hub_score"])
def paired(ax, key, title, letter, dx=-0.5):
    d = w[key][["C1D1", "C1D15"]].dropna()
    for _, r in d.iterrows():
        ax.plot([0, 1], [r.C1D1, r.C1D15], color="#C9C9C9", lw=0.6, zorder=1)
    for j, c in enumerate(["C1D1", "C1D15"]):
        ax.scatter([j] * len(d), d[c], s=14, color=[CAT[0], CAT[1]][j], zorder=3, lw=0.3, edgecolor="white")
        ax.plot([j - 0.18, j + 0.18], [d[c].median()] * 2, color=INK, lw=1.6, zorder=4)
    p = stats.wilcoxon(d.C1D1, d.C1D15).pvalue
    ax.set_xticks([0, 1]); ax.set_xticklabels(["C1D1\n(anastrozole)", "C1D15\n(+ palbociclib)"], fontsize=6)
    ax.set_xlim(-0.4, 1.4); ax.set_ylabel("Score ($z$)")
    ax.set_title(f"{title}\nn = {len(d)} paired, $P$ = {p:.3g}", fontsize=7.2)
    panel(ax, letter, dx=dx)
paired(fig.add_subplot(gs[1, 0]), "e2f_score", "Hallmark E2F targets", "C")
paired(fig.add_subplot(gs[1, 1]), "resist_score", "Resistance score (378 genes)", "D")
paired(fig.add_subplot(gs[1, 2]), "ifn_hub_score", "Interferon hub score", "E")
save(fig, "Fig6_tme_neopalana")
print(hw[hw.FDR < 0.05].head(6).round(3).to_string(index=False))
print(json.dumps(npr, indent=1)[:600])
