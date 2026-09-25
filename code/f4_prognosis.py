import numpy as np, pandas as pd, json, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fig_style import *
from surv_utils import km, logrank, time_auc
O = "/home/claude/an/out"
T = pd.read_csv(f"{O}/tcga_cohort.csv", index_col=0); M = pd.read_csv(f"{O}/metabric_cohort.csv", index_col=0)
res = json.load(open(f"{O}/survival_results.json")); coef = pd.read_csv(f"{O}/risk_coefficients.csv", index_col=0).iloc[:, 0]
cox = pd.read_csv(f"{O}/cox_uni_multi.csv"); cal = pd.read_csv(f"{O}/calibration_5y.csv"); dca = pd.read_csv(f"{O}/dca_5y.csv")
uni = pd.read_csv(f"{O}/univariable_gene_cox.csv"); lcv = pd.read_csv(f"{O}/lassocox_cv.csv")

fig = plt.figure(figsize=(7.2, 7.8))
gs = fig.add_gridspec(3, 3, hspace=0.72, wspace=0.52, height_ratios=[1, 1, 1])

# A LASSO-Cox CV + coefficients
axa = fig.add_subplot(gs[0, 0])
axa.plot(np.log10(lcv["lambda"]), lcv.cv_deviance, color=CAT[0], lw=1.3)
j = lcv.cv_deviance.idxmin(); axa.axvline(np.log10(lcv["lambda"][j]), ls="--", color=UP, lw=1)
axa.set_xlabel("log$_{10}$ $\\lambda$"); axa.set_ylabel("CV partial-likelihood deviance")
axa.set_title("LASSO-Cox in TCGA", fontsize=7.5)
ins = axa.inset_axes([0.45, 0.58, 0.52, 0.38])
cs = coef.sort_values()
ins.barh(np.arange(len(cs)), cs.values, color=[UP if v > 0 else DOWN for v in cs.values], height=0.6)
ins.set_yticks(np.arange(len(cs))); ins.set_yticklabels(cs.index, fontsize=5.2)
ins.axvline(0, color=INK2, lw=0.5); ins.tick_params(labelsize=5); ins.set_xlabel("coef.", fontsize=5.4, labelpad=1)
panel(axa, "A", dx=-0.34)

# B,C Kaplan-Meier
def kmplot(ax, D, tcol, ecol, title, xmax):
    for g, c, lbl in [(0, LOW, "Low risk"), (1, HIGH, "High risk")]:
        d = D[D.high == g].dropna(subset=[tcol, ecol])
        k = km(d[tcol].values, d[ecol].values.astype(int))
        ax.step(k.time, k.surv, where="post", color=c, label=f"{lbl} (n={len(d)})")
    d = D.dropna(subset=[tcol, ecol])
    chi, p = logrank(d[tcol].values, d[ecol].values.astype(int), d.high.values)
    ax.text(0.04, 0.08, f"log-rank $P$ = {p:.3g}", transform=ax.transAxes, fontsize=6.2)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.02); ax.set_xlabel("Years"); ax.set_ylabel("Event-free probability")
    ax.set_title(title, fontsize=7.5); ax.legend(fontsize=6, loc="lower left", bbox_to_anchor=(0.02, 0.16))
kmplot(fig.add_subplot(gs[0, 1]), T, "time", "event", "TCGA ER+/HER2$-$\nprogression-free interval", 15)
panel(fig.axes[-1], "B", dx=-0.3)
kmplot(fig.add_subplot(gs[0, 2]), M, "time", "event", "METABRIC ER+/HER2$-$\nrelapse-free survival", 25)
panel(fig.axes[-1], "C", dx=-0.3)

# D time-dependent AUC
axd = fig.add_subplot(gs[1, 0])
for (nm, D, yrs, c) in [("TCGA (PFI)", T, [3, 5, 8], CAT[0]), ("METABRIC (RFS)", M, [3, 5, 10], CAT[1])]:
    k = "TCGA_tAUC" if nm.startswith("TCGA") else "METABRIC_tAUC"
    a = [res[k][str(y)] for y in yrs]
    axd.errorbar(yrs, [x[0] for x in a], yerr=[[x[0] - x[1] for x in a], [x[2] - x[0] for x in a]],
                 fmt="o-", color=c, ms=4, capsize=2.5, lw=1.2, label=nm)
axd.axhline(0.5, color=MUTED, ls=(0, (3, 3)), lw=0.8)
axd.set_xlabel("Years"); axd.set_ylabel("Time-dependent AUC"); axd.set_ylim(0.3, 0.85)
axd.set_title("Discrimination over time", fontsize=7.5); axd.legend(fontsize=6, loc="upper left")
panel(axd, "D", dx=-0.34)

# E forest plot
axe = fig.add_subplot(gs[1, 1:])
lab = {"risk_sd": "Risk score (per SD)", "age": "Age (per year)", "stage_II": "Stage II vs I",
       "stage_III-IV": "Stage III–IV vs I", "grade_3": "Grade 3 vs 1–2"}
rows = []
for ch in ["TCGA", "METABRIC"]:
    for m in ["Univariable", "Multivariable"]:
        d = cox[(cox.cohort == ch) & (cox.model == m)]
        for _, r in d.iterrows():
            if r.term in lab: rows.append((f"{ch} — {m[:5]}.", lab[r.term], r.HR, r.lower, r.upper, r.p))
F = pd.DataFrame(rows, columns=["grp", "term", "HR", "lo", "hi", "p"])
F = F[F.term == "Risk score (per SD)"].reset_index(drop=True)
y = np.arange(len(F))[::-1]
axe.errorbar(F.HR, y, xerr=[F.HR - F.lo, F.hi - F.HR], fmt="s", color=INK, ms=4.5, capsize=3, lw=1.1)
axe.axvline(1, color=MUTED, ls=(0, (3, 3)), lw=0.8)
axe.set_yticks(y); axe.set_yticklabels(F.grp, fontsize=6.2)
for i, r in F.iterrows():
    axe.text(F.hi.max() * 1.06, y[i], f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f})   $P$={r.p:.3g}", fontsize=5.8, va="center")
axe.set_xlim(0.6, F.hi.max() * 1.9); axe.set_xlabel("Hazard ratio per SD of risk score (95% CI)")
axe.set_title("Risk score adjusted for age, stage and grade", fontsize=7.5)
panel(axe, "E", dx=-0.2)

# F calibration
axf = fig.add_subplot(gs[2, 0])
axf.plot([0, 0.45], [0, 0.45], color=MUTED, ls=(0, (3, 3)), lw=0.8)
axf.plot(cal.predicted, cal.observed, "o-", color=CAT[0], ms=4.5, mfc="white", mew=1.1)
axf.set_xlabel("Predicted 5-year relapse risk"); axf.set_ylabel("Observed (Kaplan–Meier)")
axf.set_title("Calibration, METABRIC", fontsize=7.5); axf.set_xlim(0.08, 0.42); axf.set_ylim(0.08, 0.42)
panel(axf, "F", dx=-0.34)

# G decision curve
axg = fig.add_subplot(gs[2, 1])
axg.plot(dca.threshold, dca.clinical, color=CAT[0], label="Clinical model")
axg.plot(dca.threshold, dca.full, color=CAT[1], ls="--", label="Clinical + signature")
axg.plot(dca.threshold, dca.treat_all, color=MUTED, lw=0.9, label="Treat all")
axg.axhline(0, color=INK2, lw=0.7, label="Treat none")
axg.set_xlabel("Threshold probability"); axg.set_ylabel("Net benefit")
axg.set_title("Decision curve, 5 years", fontsize=7.5); axg.legend(fontsize=5.8, loc="upper right")
axg.set_ylim(-0.03, 0.22)
panel(axg, "G", dx=-0.34)

# H per-gene univariable HRs in both cohorts
axh = fig.add_subplot(gs[2, 2])
u = uni.pivot(index="gene", columns="cohort", values="HR")
q = uni.pivot(index="gene", columns="cohort", values="p")
sig = (q < 0.05)
axh.axhline(1, color=MUTED, lw=0.7, ls=(0, (3, 3))); axh.axvline(1, color=MUTED, lw=0.7, ls=(0, (3, 3)))
axh.scatter(u.TCGA, u.METABRIC, s=18, c=[UP if (sig.TCGA[g] or sig.METABRIC[g]) else "#BBBBBB" for g in u.index],
            lw=0.4, edgecolor="white")
for g in u.index:
    if sig.TCGA[g] or sig.METABRIC[g]:
        axh.annotate(g, (u.TCGA[g], u.METABRIC[g]), fontsize=5.4, xytext=(3, 2), textcoords="offset points")
axh.set_xlabel("HR per SD, TCGA (PFI)"); axh.set_ylabel("HR per SD, METABRIC (RFS)")
axh.set_title("Per-gene prognostic effect", fontsize=7.5)
panel(axh, "H", dx=-0.34)
save(fig, "Fig5_prognosis")
print(json.dumps({k: v for k, v in res.items() if k.endswith(("primary", "OS"))}, indent=1, default=float)[:1200])
