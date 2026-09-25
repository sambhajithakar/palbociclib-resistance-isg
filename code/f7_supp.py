import numpy as np, pandas as pd, json, matplotlib.pyplot as plt
from fig_style import *
from surv_utils import km, logrank
O = "/home/claude/an/out"
core = json.load(open(f"{O}/ml_selection.json"))["core"]
M = pd.read_pickle(f"{O}/merged_centered.pkl"); meta = pd.read_csv(f"{O}/meta.csv").set_index("sample").loc[M.columns]

# --- S1: core signature expression ---
fig, ax = plt.subplots(figsize=(7.2, 3.4))
order = sorted(core)
X = M.loc[order]; X = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
o = np.lexsort((meta.dataset.values, (meta.group != "Resistant").values))
im = ax.imshow(X.values[:, o], cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
nres = (meta.group == "Resistant").sum()
ax.axvline(nres - 0.5, color="black", lw=1.2)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=6)
ax.set_xticks([nres / 2, nres + (len(o) - nres) / 2]); ax.set_xticklabels(["Resistant (n=49)", "Sensitive (n=25)"], fontsize=7)
ax.set_title("Core signature expression across discovery samples", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.015); cb.set_label("$z$-score", fontsize=6.5); cb.ax.tick_params(labelsize=6)
for s in ax.spines.values(): s.set_visible(True)
save(fig, "FigS1_core_expression")

# --- S2: overall survival ---
T = pd.read_csv(f"{O}/tcga_cohort.csv", index_col=0); Mb = pd.read_csv(f"{O}/metabric_cohort.csv", index_col=0)
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
for ax, (D, ttl, xmax) in zip(axes, [(T, "TCGA ER+/HER2$-$ — overall survival", 15),
                                     (Mb, "METABRIC ER+/HER2$-$ — overall survival", 25)]):
    for g, c, lbl in [(0, LOW, "Low risk"), (1, HIGH, "High risk")]:
        d = D[D.high == g].dropna(subset=["os_time", "os"])
        k = km(d.os_time.values, d.os.values.astype(int))
        ax.step(k.time, k.surv, where="post", color=c, label=f"{lbl} (n={len(d)})")
    d = D.dropna(subset=["os_time", "os"])
    chi, p = logrank(d.os_time.values, d.os.values.astype(int), d.high.values)
    ax.text(0.04, 0.08, f"log-rank $P$ = {p:.3g}", transform=ax.transAxes, fontsize=7)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.02); ax.set_xlabel("Years"); ax.set_ylabel("Overall survival")
    ax.set_title(ttl, fontsize=8); ax.legend(fontsize=7, loc="lower left", bbox_to_anchor=(0.02, 0.16))
fig.tight_layout()
save(fig, "FigS2_overall_survival")
