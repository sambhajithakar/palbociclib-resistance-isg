import numpy as np, pandas as pd, matplotlib.pyplot as plt, pickle
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa
from fig_style import *
D = "/home/claude/an/dock"
s1 = pd.read_csv(f"{D}/stage1_scores.csv").dropna(subset=["vina_kcal_mol"])
s2 = pd.read_csv(f"{D}/stage2_top.csv")
s2 = s2.drop_duplicates("smiles", keep="first") if "smiles" in s2 else s2
props = pd.read_csv(f"{D}/library_properties.csv")
val = pd.read_csv(f"{D}/redock_validation.csv")
ref = np.load(f"{D}/ref_ptyr.npy"); P = np.load(f"{D}/phosphate_center.npy")
CTRL = "phosphotyrosine (crystal ligand)"

fig = plt.figure(figsize=(7.2, 7.2))
gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.42, height_ratios=[1, 1, 1.1])

# A redocking validation
axa = fig.add_subplot(gs[0, 0])
x = np.arange(len(val))
b = axa.bar(x, val.mean_nn_dist, color=[CAT[2] if v < 2 else "#C4C4C4" for v in val.mean_nn_dist], width=0.5)
axa.axhline(2.0, color=UP, ls="--", lw=1.1)
axa.text(len(val) - 0.5, 2.1, "2 Å acceptance", ha="right", fontsize=6, color=UP)
for r, v in zip(b, val.mean_nn_dist): axa.text(r.get_x() + r.get_width() / 2, v + 0.08, f"{v:.2f}", ha="center", fontsize=6.2)
axa.set_xticks(x); axa.set_xticklabels([f"{int(w)} Å\n{c}" for w, c in zip(val.box, val.centre)], fontsize=6)
axa.set_ylabel("Mean distance to crystal\nphosphotyrosine atoms (Å)")
axa.set_title("Redocking validation", fontsize=7.5); axa.set_ylim(0, max(val.mean_nn_dist) * 1.3)
panel(axa, "A", dx=-0.34)

# B score distribution
axb = fig.add_subplot(gs[0, 1])
axb.hist(s1.vina_kcal_mol, bins=40, color="#BFD4E8", edgecolor="white", linewidth=0.3)
for nm, c, dy in [(CTRL, CAT[2], 0), ("palbociclib", CAT[3], 0)]:
    v = s1.loc[s1.compound == nm, "vina_kcal_mol"]
    if len(v): axb.axvline(v.values[0], color=c, lw=1.4, ls="--")
best = s2[~s2.compound.isin([CTRL, "palbociclib"])].vina_mean.min()
axb.set_xlabel("Stage-1 Vina score (kcal/mol)"); axb.set_ylabel("Approved drugs")
axb.set_title(f"Stage-1 screen of {len(s1)-2:,} approved drugs", fontsize=7.5)
axb.legend(handles=[Line2D([], [], color=CAT[2], ls="--", label="phosphotyrosine"),
                    Line2D([], [], color=CAT[3], ls="--", label="palbociclib")], fontsize=6, loc="upper left")
panel(axb, "B", dx=-0.3)

# C top compounds
axc = fig.add_subplot(gs[1, :])
t = s2.sort_values("vina_mean").head(16).iloc[::-1]
y = np.arange(len(t))
cols = [CAT[2] if c == CTRL else CAT[3] if c == "palbociclib" else CAT[0] for c in t.compound]
axc.barh(y, -t.vina_mean, xerr=t.vina_sd, color=cols, height=0.62,
         error_kw=dict(lw=0.8, capsize=2, ecolor=INK2))
axc.set_yticks(y); axc.set_yticklabels([c[:34] for c in t.compound], fontsize=6)
axc.set_xlabel("$-$Vina score (kcal/mol), mean ± SD of three runs")
axc.set_title("Stage-2 re-docking at exhaustiveness 32", fontsize=7.5)
axc.set_xlim(0, (-t.vina_mean + t.vina_sd).max() * 1.1)
axc.legend(handles=[Line2D([], [], color=CAT[0], lw=4, label="Approved drug"),
                    Line2D([], [], color=CAT[2], lw=4, label="Native ligand"),
                    Line2D([], [], color=CAT[3], lw=4, label="Palbociclib")], fontsize=6, loc="lower right")
panel(axc, "C", dx=-0.16)

# D contact residues of the leading compounds
axd = fig.add_subplot(gs[2, 0])
lead = s2.sort_values("vina_mean").head(8)
allres = []
for c in lead.contacts.fillna(""):
    allres += [r for r in c.split(";") if r]
cnt = pd.Series(allres).value_counts().head(12).iloc[::-1]
axd.barh(np.arange(len(cnt)), cnt.values, color=CAT[0], height=0.64)
axd.set_yticks(np.arange(len(cnt))); axd.set_yticklabels(cnt.index, fontsize=6)
axd.set_xlabel("Leading compounds in contact (of 8)")
axd.set_title("Contacted residues (< 4 Å)", fontsize=7.5)
axd.set_xlim(0, 8.6)
panel(axd, "D", dx=-0.34)

# E property space
axe = fig.add_subplot(gs[2, 1])
scr = props[props.name.isin(s1.compound)]
axe.scatter(scr.MW, scr.LogP, s=5, c="#DCDCDC", lw=0, label="Screened library")
m = props.set_index("name")
sel = [c for c in s2.compound if c in m.index]
axe.scatter(m.loc[sel].MW, m.loc[sel].LogP, s=26, c=CAT[0], lw=0.4, edgecolor="white", label="Stage-2 compounds")
if CTRL in m.index:
    axe.scatter(m.loc[[CTRL]].MW, m.loc[[CTRL]].LogP, s=42, c=CAT[2], marker="D", lw=0.4, edgecolor="white", label="phosphotyrosine")
axe.axhline(5, color=MUTED, ls=(0, (3, 3)), lw=0.8)
axe.set_xlabel("Molecular weight"); axe.set_ylabel("cLogP")
axe.set_title("Property space", fontsize=7.5); axe.legend(fontsize=6, loc="lower right")
panel(axe, "E", dx=-0.3)
save(fig, "Fig7_docking")
print(s2.head(12)[["compound", "vina_mean", "vina_sd", "n_contacts", "hbond_residues"]].to_string(index=False))
