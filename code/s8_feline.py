"""External validation in the FELINE trial (GSE158724).

ER+ patients randomised to neoadjuvant letrozole with placebo (arm A) or with
ribociclib (arms B, C), biopsied at day 0 (start), day 14 (mid) and day 180 (end),
with a clinical response annotation. Cancer cells are pseudobulked per biopsy and
the resistance programme, the interferon hub set, Hallmark E2F targets and the
four-gene risk score are applied.

Usage:  python s8_feline.py [10x|icell8]
"""
import gzip, json, sys, numpy as np, pandas as pd
from scipy import stats

R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"
O = "/home/claude/an/out"
PLAT = sys.argv[1] if len(sys.argv) > 1 else "10x"
MIN_CELLS = 20

meta = pd.read_csv(f"{R}/FELINE_sample_metadata.tsv", sep="\t")
ct = pd.read_csv(f"{R}/FELINE_{PLAT}_celltypes.txt.gz", sep="\t")
ct.columns = ["cell", "type"]
cancer = set(ct.cell[ct.type.str.contains("Cancer", case=False, na=False)])
print(f"{PLAT}: {len(ct)} annotated cells, {len(cancer)} cancer cells")

# ---------------------------------------------------------------- gene sets
rob = pd.read_csv(f"{O}/robust_DEGs.csv")
up = list(rob.gene[rob.direction == "Up"])
dn = list(rob.gene[rob.direction == "Down"])
hub = ["STAT1", "ISG15", "IRF1", "IFIT3", "IRF7", "IFI44", "DDX60", "IFIT1",
       "RSAD2", "IFIT2", "OASL", "IFI6", "HELZ2", "CMPK2"]
coef = pd.read_csv(f"{O}/risk_coefficients.csv", index_col=0).iloc[:, 0]
E2F = [l for l in open(f"{R}/MSigDB_Hallmark_2020.gmt") if l.startswith("E2F Targets")][0]
E2F = [g.strip() for g in E2F.split("\t")[2:] if g.strip()]
need = set(up) | set(dn) | set(hub) | set(coef.index) | set(E2F)

# ------------------------------------------------- stream counts, keep needed
import os
counts = f"/home/claude/FELINE_{PLAT}_counts.txt.gz"
if not os.path.exists(counts):
    counts = f"{R}/FELINE_{PLAT}_counts.txt.gz"

# read in row blocks with the C parser: accumulate the library size of every cell
# and retain only the gene rows the scores need
with gzip.open(counts, "rt") as _f:
    _hdr = _f.readline().rstrip("\n").split("\t")
dt = {c: np.float32 for c in _hdr[1:]}
dt[_hdr[0]] = str

kept, tot, rows, n = None, None, {}, 0
for blk in pd.read_csv(counts, sep="\t", index_col=0, chunksize=400,
                       dtype=dt, engine="c"):
    if kept is None:
        kept = [c for c in blk.columns]
        keep_mask = np.array([c in cancer for c in kept])
        kept = [c for c, k in zip(kept, keep_mask) if k]
        tot = np.zeros(len(kept), dtype=np.float64)
    V = blk.values[:, keep_mask]
    tot += V.sum(axis=0)
    n += V.shape[0]
    hit = [i for i, g in enumerate(blk.index) if g in need]
    for i in hit:
        rows[blk.index[i]] = np.array(V[i], copy=True)
    if n % 2000 == 0:
        print(f"  {n} genes read, {len(rows)} retained", flush=True)
print(f"genes read: {n}, retained: {len(rows)}; cancer cells in matrix: {len(kept)}")

X = pd.DataFrame(rows, index=kept).T
cpm = np.log2(X.div(pd.Series(tot, index=kept), axis=1) * 1e6 + 1)

# -------------------------------------------- pseudobulk per biopsy (patient_tp)
samp = pd.Series(["_".join(c.split("_")[:2]) for c in kept], index=kept)
ncell = samp.value_counts()
pb = cpm.T.groupby(samp).mean().T
pb = pb[[s for s in pb.columns if ncell[s] >= MIN_CELLS]]
print(f"biopsies with >={MIN_CELLS} cancer cells: {pb.shape[1]}")

Z = pb.sub(pb.mean(1), axis=0).div(pb.std(1) + 1e-9, axis=0)


def score(genes):
    g = [x for x in genes if x in Z.index]
    return (Z.loc[g].mean() if g else pd.Series(np.nan, index=Z.columns)), len(g)


s_up, n_up = score(up)
s_dn, n_dn = score(dn)
s_hub, n_hub = score(hub)
s_e2f, n_e2f = score(E2F)
res = pd.DataFrame({"resist_score": s_up - s_dn, "ifn_hub": s_hub, "e2f": s_e2f})
cg = [g for g in coef.index if g in Z.index]
res["risk4"] = (Z.loc[cg].T @ coef[cg]) if cg else np.nan
res["n_cells"] = ncell.reindex(res.index).values
res = res.join(meta.drop_duplicates("title").set_index("title")[
    ["patient", "arm", "treatment", "response", "timepoint", "days", "schedule"]], how="left")
res["ribo"] = res.arm.isin(["B", "C"])
res = res.dropna(subset=["patient"])
res.to_csv(f"{O}/feline_{PLAT}_scores.csv")
print(f"genes matched: up {n_up}/{len(up)}, down {n_dn}/{len(dn)}, "
      f"hub {n_hub}/{len(hub)}, E2F {n_e2f}/{len(E2F)}, risk {len(cg)}/{len(coef)}")
print(res.groupby(["arm", "timepoint"]).size().to_dict())

SCORES = ["resist_score", "ifn_hub", "e2f", "risk4"]
out = {}


def wilcox_paired(sub, sc, a, b, label):
    w = sub.pivot_table(index="patient", columns="timepoint", values=sc)
    if a not in w.columns or b not in w.columns:
        return
    x = w[[a, b]].dropna()
    if len(x) < 5:
        return
    st = stats.wilcoxon(x[a], x[b])
    out[label] = dict(n=int(len(x)), median_change=float((x[b] - x[a]).median()),
                      p=float(st.pvalue))


# (1) on-treatment change, ribociclib vs letrozole alone
for lbl, sub in [("ribociclib", res[res.ribo]), ("letrozole", res[~res.ribo])]:
    for a, b in [("start", "mid"), ("start", "end")]:
        for sc in SCORES:
            wilcox_paired(sub, sc, a, b, f"change|{lbl}|{sc}|{a}->{b}")

# (2) baseline score vs response, by arm
for lbl, sub in [("ribociclib", res[res.ribo]), ("letrozole", res[~res.ribo])]:
    base = sub[(sub.timepoint == "start")].dropna(subset=["response"])
    for sc in SCORES:
        r = base[base.response == "Responder"][sc].dropna()
        nr = base[base.response == "Non-responder"][sc].dropna()
        if len(r) >= 4 and len(nr) >= 4:
            u = stats.mannwhitneyu(nr, r, alternative="two-sided")
            out[f"baseline|{lbl}|{sc}"] = dict(
                n_resp=int(len(r)), n_non=int(len(nr)),
                median_resp=float(r.median()), median_non=float(nr.median()),
                AUC_high_predicts_nonresponse=float(u.statistic / (len(r) * len(nr))),
                p=float(u.pvalue))

# (3) day-180 score vs response, ribociclib arms
end = res[(res.timepoint == "end") & res.ribo].dropna(subset=["response"])
for sc in SCORES:
    r = end[end.response == "Responder"][sc].dropna()
    nr = end[end.response == "Non-responder"][sc].dropna()
    if len(r) >= 4 and len(nr) >= 4:
        u = stats.mannwhitneyu(nr, r, alternative="two-sided")
        out[f"day180|ribociclib|{sc}"] = dict(
            n_resp=int(len(r)), n_non=int(len(nr)),
            AUC_high_predicts_nonresponse=float(u.statistic / (len(r) * len(nr))),
            p=float(u.pvalue))

# (4) on-treatment change vs response (ribociclib): does the programme rise in
#     patients who fail?
for a, b in [("start", "mid"), ("start", "end")]:
    for sc in SCORES:
        w = res[res.ribo].pivot_table(index="patient", columns="timepoint", values=sc)
        if a not in w.columns or b not in w.columns:
            continue
        d = (w[b] - w[a]).dropna()
        rsp = res.drop_duplicates("patient").set_index("patient").response
        r = d[rsp.reindex(d.index) == "Responder"]
        nr = d[rsp.reindex(d.index) == "Non-responder"]
        if len(r) >= 4 and len(nr) >= 4:
            u = stats.mannwhitneyu(nr, r, alternative="two-sided")
            out[f"delta|ribociclib|{sc}|{a}->{b}"] = dict(
                n_resp=int(len(r)), n_non=int(len(nr)),
                median_change_resp=float(r.median()), median_change_non=float(nr.median()),
                p=float(u.pvalue))

json.dump(out, open(f"{O}/feline_{PLAT}_results.json", "w"), indent=1)
for k, v in out.items():
    print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items()})
