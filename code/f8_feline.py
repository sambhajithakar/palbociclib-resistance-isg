"""Figure S3: external validation of the resistance programme in the FELINE trial."""
import json, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
from fig_style import CAT, INK, INK2, MUTED, panel, save

O = "/home/claude/an/out"
res = pd.read_csv(f"{O}/feline_10x_scores.csv", index_col=0)
R = json.load(open(f"{O}/feline_10x_results.json"))
res["arm_lbl"] = np.where(res.ribo, "Letrozole + ribociclib", "Letrozole")
TP = {"start": "Day 0", "mid": "Day 14", "end": "Day 180"}

fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.42)


def strip(ax, sub, sc, ylab, title):
    grps = ["Responder", "Non-responder"]
    for i, g in enumerate(grps):
        v = sub[sub.response == g][sc].dropna().values
        x = np.random.default_rng(3).normal(i, 0.055, len(v))
        ax.scatter(x, v, s=16, color=CAT[i], alpha=.85, linewidth=0, zorder=3)
        ax.hlines(np.median(v), i - .22, i + .22, color=INK, lw=1.6, zorder=4)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Responder", "Non-\nresponder"])
    ax.set_xlim(-.5, 1.5); ax.set_ylabel(ylab); ax.set_title(title, pad=9)
    ax.grid(axis="y", zorder=0)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * .20)


# A-C: baseline scores by response, ribociclib arms
base = res[(res.timepoint == "start") & res.ribo]
for k, (sc, lab) in enumerate([("resist_score", "Resistance score"),
                               ("ifn_hub", "Interferon hub score"),
                               ("risk4", "Four-gene risk score")]):
    ax = fig.add_subplot(gs[0, k])
    strip(ax, base, sc, lab, "Baseline, ribociclib arms")
    d = R.get(f"baseline|ribociclib|{sc}", {})
    if d:
        ax.text(.5, .985, f"AUC {d['AUC_high_predicts_nonresponse']:.2f}, P = {d['p']:.2f}",
                transform=ax.transAxes, ha="center", va="top", fontsize=6.5, color=INK2)
    panel(ax, "ABC"[k])

# D: ROC-free AUC summary across scores and settings
ax = fig.add_subplot(gs[1, 0])
keys = [("baseline|ribociclib|resist_score", "Resistance, day 0"),
        ("baseline|ribociclib|ifn_hub", "IFN hub, day 0"),
        ("baseline|ribociclib|risk4", "Risk score, day 0"),
        ("day180|ribociclib|resist_score", "Resistance, day 180"),
        ("day180|ribociclib|ifn_hub", "IFN hub, day 180"),
        ("day180|ribociclib|risk4", "Risk score, day 180")]
lab = [l for k_, l in keys if k_ in R]
val = [R[k_]["AUC_high_predicts_nonresponse"] for k_, _ in keys if k_ in R]
y = np.arange(len(val))
ax.barh(y, val, color=CAT[0], height=.62)
ax.axvline(.5, color=MUTED, lw=.9, ls="--")
ax.set_yticks(y); ax.set_yticklabels(lab); ax.invert_yaxis()
ax.set_xlim(0, 1); ax.set_xlabel("AUC for non-response")
ax.set_title("Discrimination of clinical response", pad=9)
ax.grid(axis="x"); panel(ax, "D")

# E-F: on-treatment change, day 0 to day 14, by arm
for k, (sc, lab2) in enumerate([("e2f", "Hallmark E2F targets"),
                                ("resist_score", "Resistance score")]):
    ax = fig.add_subplot(gs[1, 1 + k]); ptxt = []
    for i, (arm, sub) in enumerate([("Letrozole", res[~res.ribo]),
                                    ("Letrozole\n+ ribociclib", res[res.ribo])]):
        w = sub.pivot_table(index="patient", columns="timepoint", values=sc)
        if "start" not in w or "mid" not in w:
            continue
        x = w[["start", "mid"]].dropna()
        for _, r in x.iterrows():
            ax.plot([i - .16, i + .16], [r["start"], r["mid"]],
                    color=CAT[i], lw=.8, alpha=.55, zorder=2)
        ax.scatter([i - .16] * len(x), x["start"], s=13, color=CAT[i], zorder=3, linewidth=0)
        ax.scatter([i + .16] * len(x), x["mid"], s=13, color=CAT[i], zorder=3,
                   facecolor="white", edgecolor=CAT[i], linewidth=.9)
        key = f"change|{'ribociclib' if i else 'letrozole'}|{sc}|start->mid"
        if key in R:
            ptxt.append((i, f"P = {R[key]['p']:.3f}"))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Letrozole", "Letrozole\n+ ribociclib"])
    ax.set_xlim(-.5, 1.5); ax.set_ylabel(lab2)
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * .16)
    ax.set_title("Day 0 (filled) to day 14 (open)", pad=9)
    for xi, t in ptxt:
        ax.text(xi, .99, t, ha="center", va="top", transform=ax.get_xaxis_transform(),
                fontsize=6.3, color=INK2)
    ax.grid(axis="y"); panel(ax, "EF"[k])

save(fig, "figS3")
