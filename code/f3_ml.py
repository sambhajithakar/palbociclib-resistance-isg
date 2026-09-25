import numpy as np, pandas as pd, json, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, roc_auc_score
from fig_style import *
O = "/home/claude/an/out"
sel = json.load(open(f"{O}/ml_selection.json")); core = sel["core"]; sets = sel["sets"]
lcv = pd.read_csv(f"{O}/lasso_cv.csv"); lpath = pd.read_csv(f"{O}/lasso_path.csv", index_col=0)
scv = pd.read_csv(f"{O}/svmrfe_cv.csv"); rf = pd.read_csv(f"{O}/rf_boruta.csv")
auc = pd.read_csv(f"{O}/ml_auc.csv"); pred = pd.read_csv(f"{O}/ml_predictions.csv")
gbm = pd.read_csv(f"{O}/gbm_importance.csv"); gpred = pd.read_csv(f"{O}/gbm_predictions.csv")
M = pd.read_pickle(f"{O}/merged_centered.pkl"); meta = pd.read_csv(f"{O}/meta.csv").set_index("sample")

fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(3, 2, hspace=0.68, wspace=0.42, height_ratios=[1, 1, 1.05])

# A LASSO
axa = fig.add_subplot(gs[0, 0])
for g in lpath.columns:
    axa.plot(np.log10(lpath.index), lpath[g], lw=0.7, color=UP if g in core else "#D6D6D6", zorder=2 if g in core else 1)
axa.set_xlabel("log$_{10}$ $C$ (inverse penalty)"); axa.set_ylabel("Coefficient")
axa.set_title("LASSO logistic regression", fontsize=7.5)
ax2 = axa.twinx(); ax2.plot(np.log10(lcv.C), lcv.cv_logloss, color=CAT[0], lw=1.4, ls="--")
ax2.set_ylabel("10-fold CV log-loss", color=CAT[0], fontsize=7); ax2.tick_params(labelsize=6, colors=CAT[0])
ax2.spines["right"].set_visible(True); ax2.spines["right"].set_color(CAT[0]); ax2.grid(False)
axa.legend(handles=[Line2D([], [], color=UP, lw=1.2, label="Selected"),
                    Line2D([], [], color="#D6D6D6", lw=1.2, label="Not selected"),
                    Line2D([], [], color=CAT[0], lw=1.2, ls="--", label="CV log-loss")],
           fontsize=5.6, loc="upper left")
panel(axa, "A")

# B SVM-RFE
axb = fig.add_subplot(gs[0, 1])
axb.plot(scv.k, scv.cv_acc, "o-", color=CAT[0], ms=3.4, mfc="white", mew=0.9)
best = int(scv.k[scv.cv_acc.idxmax()])
axb.axvline(best, ls=(0, (3, 3)), color=UP, lw=0.9)
axb.annotate(f"{best} genes\n(accuracy {scv.cv_acc.max():.2f})", (best, scv.cv_acc.max()),
             xytext=(10, -14), textcoords="offset points", fontsize=6, color=UP)
axb.set_xlabel("Number of genes retained"); axb.set_ylabel("10-fold CV accuracy")
axb.set_title("SVM recursive feature elimination", fontsize=7.5)
panel(axb, "B")

# C Boruta / RF
axc = fig.add_subplot(gs[1, 0])
r = rf.nlargest(18, "gini").iloc[::-1]
cols = [UP if p < 0.05 else "#C9C9C9" for p in r.boruta_padj]
axc.barh(np.arange(len(r)), r.gini, color=cols, height=0.64)
axc.set_yticks(np.arange(len(r))); axc.set_yticklabels(r.gene, fontsize=5.8)
axc.set_xlabel("Random-forest importance (Gini)"); axc.set_title("Random forest with Boruta", fontsize=7.5)
axc.legend(handles=[Line2D([], [], color=UP, lw=4, label="Boruta-confirmed"),
                    Line2D([], [], color="#C9C9C9", lw=4, label="Rejected")], fontsize=5.6, loc="lower right")
panel(axc, "C", dx=-0.38)

# D selection matrix
axd = fig.add_subplot(gs[1, 1])
meths = ["LASSO", "SVM-RFE", "RF-Boruta"]
gtop = list(gbm.nlargest(15, "importance").gene)
order = sorted(core, key=lambda g: (-sum(g in sets[m] for m in meths), g))
mat = np.array([[1 if g in sets[m] else 0 for g in order] for m in meths] + [[1 if g in gtop else 0 for g in order]])
axd.imshow(mat, cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
for i in range(mat.shape[0] + 1): axd.axhline(i - 0.5, color="white", lw=1.2)
for j in range(mat.shape[1] + 1): axd.axvline(j - 0.5, color="white", lw=1.2)
axd.set_yticks(range(4)); axd.set_yticklabels(meths + ["GBM top-15"], fontsize=6)
axd.set_xticks(range(len(order))); axd.set_xticklabels(order, rotation=90, fontsize=5.6)
axd.set_title(f"Core signature: {len(core)} genes selected by ≥ 2 of 3 algorithms", fontsize=7.2)
for s in axd.spines.values(): s.set_visible(False)
panel(axd, "D", dx=-0.16, dy=1.12)

# E ROC
axe = fig.add_subplot(gs[2, 0])
y = pred.y.values
for lbl, p, c, ls in [("Combined, 10-fold CV", pred.p_cv.values, CAT[0], "-"),
                      ("Combined, leave-one-dataset-out", pred.p_lodo.values, CAT[1], "--"),
                      ("Gradient boosting, CV", gpred.gbm_core.values, CAT[2], "-.")]:
    fpr, tpr, _ = roc_curve(y, p)
    axe.plot(fpr, tpr, color=c, ls=ls, label=f"{lbl} (AUC {roc_auc_score(y, p):.3f})")
axe.plot([0, 1], [0, 1], color=MUTED, lw=0.7, ls=(0, (3, 3)))
axe.set_xlabel("1 $-$ specificity"); axe.set_ylabel("Sensitivity")
axe.set_title("Discrimination in the discovery set", fontsize=7.5)
axe.legend(fontsize=5.6, loc="lower right"); axe.set_xlim(-0.02, 1.02); axe.set_ylim(-0.02, 1.02)
panel(axe, "E")

# F nested versus apparent validation
axf = fig.add_subplot(gs[2, 1])
nst = pd.read_csv(f"{O}/nested_lodo.csv")
x = np.arange(len(nst))
b = axf.bar(x, nst.AUC, color=CAT[1], width=0.55)
for r, v, n in zip(b, nst.AUC, nst.n_test):
    axf.text(r.get_x() + r.get_width() / 2, v + 0.02, f"{v:.2f}", ha="center", fontsize=6.2)
    axf.text(r.get_x() + r.get_width() / 2, 0.03, f"n={n}", ha="center", fontsize=5.6, color="white")
axf.axhline(1.0, color=UP, lw=1.2, ls="--")
axf.text(len(nst) - 0.5, 1.02, "apparent AUC (non-nested) = 1.00", ha="right", fontsize=5.8, color=UP)
axf.axhline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)))
axf.axhline(nst.AUC.mean(), color=CAT[0], lw=1.0)
axf.text(-0.45, nst.AUC.mean() + 0.02, f"mean {nst.AUC.mean():.2f}", fontsize=5.8, color=CAT[0])
axf.set_xticks(x); axf.set_xticklabels([d.replace("GSE", "GSE ") for d in nst.held_out], fontsize=6)
axf.set_xlabel("Held-out dataset"); axf.set_ylabel("AUC in held-out dataset")
axf.set_ylim(0, 1.16); axf.set_title("Fully nested leave-one-dataset-out", fontsize=7.5)
panel(axf, "F", dx=-0.2)
save(fig, "Fig4_machine_learning")
print(auc.round(3).to_string(index=False))
