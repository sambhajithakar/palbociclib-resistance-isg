import numpy as np, pandas as pd, re, json
from stats_utils import lm_ebayes, bh, logcpm, rra
R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"
OUT = "/home/claude/an/out"; import os; os.makedirs(OUT, exist_ok=True)
LFC, FDR = 1.0, 0.05

annot = pd.read_csv(f"{R}/Human.GRCh38.p13.annot.tsv.gz", sep="\t", usecols=["GeneID", "Symbol", "GeneType"])
pc = annot[annot.GeneType == "protein-coding"].set_index("GeneID").Symbol

def ncbi(gse):
    m = pd.read_csv(f"{R}/{gse}_raw_counts.tsv.gz", sep="\t", index_col=0)
    m = m.loc[m.index.isin(pc.index)]
    m.index = pc.loc[m.index].values
    return m.groupby(level=0).sum()

# ---- sample groups (verified from GEO sample records) ----
g1 = ncbi("GSE130437")
s1 = pd.DataFrame({"sample": ["GSM3738651","GSM3738652","GSM3738653","GSM3738654","GSM3738655","GSM3738656"],
                   "group": ["Resistant"]*3 + ["Sensitive"]*3, "cell": "MCF7"})
g2 = ncbi("GSE222367")
t2 = {}
ids = [f"GSM69216{n}" for n in range(37, 79)]
# 37-39 MCF7 parental, 40-51 MCF7 palbo-R (1.2/2.4/3.6/4.8 uM), 52-54 MCF7 parental, 55-63 abema-R (excluded),
# 64-66 T47D parental, 67-78 T47D palbo-R
rows = []
for n in range(37, 79):
    gsm = f"GSM69216{n}"
    if n in (37,38,39,52,53,54): rows.append((gsm, "Sensitive", "MCF7"))
    elif 40 <= n <= 51: rows.append((gsm, "Resistant", "MCF7"))
    elif 64 <= n <= 66: rows.append((gsm, "Sensitive", "T47D"))
    elif 67 <= n <= 78: rows.append((gsm, "Resistant", "T47D"))
s2 = pd.DataFrame(rows, columns=["sample", "group", "cell"])
g3r = pd.read_csv(f"{R}/GSE229235_counts.csv.gz")
g3r = g3r[g3r.Gene.isin(set(pc.values))]
g3 = g3r.drop(columns="ID").groupby("Gene").sum()
cols3 = [c for c in g3.columns if c.startswith("ET-R") or c.startswith("PR-")]
s3 = pd.DataFrame({"sample": cols3, "group": ["Sensitive" if c.startswith("ET-R") else "Resistant" for c in cols3],
                   "cell": ["PDX-" + c.split()[0] for c in cols3]})
sets = {"GSE130437": (g1, s1), "GSE222367": (g2, s2), "GSE229235": (g3, s3)}

res_all, expr_all, meta_all = {}, {}, []
for gse, (cts, meta) in sets.items():
    miss = [x for x in meta["sample"] if x not in cts.columns]
    if miss: print(gse, "not in NCBI count matrix:", miss)
    meta = meta[meta["sample"].isin(cts.columns)].reset_index(drop=True); sets[gse] = (cts, meta)
    cts = cts[meta["sample"]]
    grp = (meta.group == "Resistant").astype(int).values
    keep = (cts >= 10).sum(1) >= min(np.bincount(grp))
    cts = cts[keep]
    y, sf = logcpm(cts)
    X = [np.ones(len(grp)), grp]
    if gse == "GSE222367":
        X.append((meta.cell == "T47D").astype(int).values)
    X = np.column_stack(X)
    lfc, t, p, d0, s02 = lm_ebayes(y.values, X, 1)
    r = pd.DataFrame({"gene": y.index, "log2FC": lfc, "t": t, "P": p, "FDR": bh(p), "AveExpr": y.mean(1).values})
    r["sig"] = np.where((r.FDR < FDR) & (r.log2FC > LFC), "Up", np.where((r.FDR < FDR) & (r.log2FC < -LFC), "Down", "NS"))
    r = r.sort_values("P")
    r.to_csv(f"{OUT}/DEG_{gse}.csv", index=False)
    res_all[gse] = r; expr_all[gse] = y
    m = meta.copy(); m["dataset"] = gse; meta_all.append(m)
    print(gse, "n=", len(grp), f"(R={grp.sum()}, S={len(grp)-grp.sum()})", "genes=", len(r),
          "Up=", (r.sig == "Up").sum(), "Down=", (r.sig == "Down").sum(), "d0=%.2f" % d0)

universe = sorted(set.intersection(*[set(r.gene) for r in res_all.values()]))
N = len(universe)
def ranked(direction):
    L = []
    for r in res_all.values():
        r = r[r.gene.isin(universe)]
        r = r[r.log2FC > 0] if direction == "up" else r[r.log2FC < 0]
        L.append(list(r.sort_values("P").gene))
    return L
ru = rra(ranked("up"), N); ru["direction"] = "Up"
rd = rra(ranked("down"), N); rd["direction"] = "Down"
ra = pd.concat([ru, rd])
lfcm = pd.DataFrame({k: v.set_index("gene").log2FC.reindex(universe) for k, v in res_all.items()})
ra["meanLFC"] = lfcm.mean(1).reindex(ra.gene).values
for k in res_all: ra[f"LFC_{k}"] = lfcm[k].reindex(ra.gene).values
ra["consistent"] = np.sign(lfcm).reindex(ra.gene).apply(lambda s: abs(s.sum()) == len(s), axis=1).values
robust = ra[(ra.score < 0.05) & (ra.meanLFC.abs() > LFC) & ra.consistent &
            (((ra.direction == "Up") & (ra.meanLFC > 0)) | ((ra.direction == "Down") & (ra.meanLFC < 0)))]
robust = robust.sort_values("score")
ra.to_csv(f"{OUT}/RRA_all.csv", index=False); robust.to_csv(f"{OUT}/robust_DEGs.csv", index=False)
print("universe", N, "robust", len(robust), "Up", (robust.direction == "Up").sum(), "Down", (robust.direction == "Down").sum())
print(robust.head(30)[["gene", "direction", "score", "meanLFC"]].to_string())

# merged expression for ML: per-dataset (per cell line/PDX line) centring = batch correction on log scale
meta = pd.concat(meta_all).reset_index(drop=True)
Z = []
for gse, y in expr_all.items():
    y = y.reindex(universe)
    m = meta[meta.dataset == gse]
    for b, mm in (m.groupby("cell") if gse != "GSE229235" else [("all", m)]):
        sub = y[mm["sample"]]
        Z.append(sub.sub(sub.mean(1), axis=0))
merged = pd.concat(Z, axis=1).dropna()
merged.to_pickle(f"{OUT}/merged_centered.pkl"); meta.to_csv(f"{OUT}/meta.csv", index=False)
pd.to_pickle(expr_all, f"{OUT}/expr_all.pkl")
print("merged", merged.shape)
