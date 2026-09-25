import gzip, re, numpy as np, pandas as pd, json
from scipy import stats
R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"; OUT = "/home/claude/an/out"
lines = gzip.open(f"{R}/GSE93204_series_matrix.txt.gz", "rt").read().split("\n")
titles = [x.strip('"') for x in next(l for l in lines if l.startswith("!Sample_title")).split("\t")[1:]]
gsm = [x.strip('"') for x in next(l for l in lines if l.startswith("!Sample_geo_accession")).split("\t")[1:]]
s = lines.index("!series_matrix_table_begin"); e = lines.index("!series_matrix_table_end")
from io import StringIO
mat = pd.read_csv(StringIO("\n".join(lines[s + 1:e])), sep="\t", index_col=0)
print(mat.shape, mat.iloc[:3, :4]); print(mat.describe().T.head(3))
gpl = [l for l in open(f"{R}/GPL6480.txt") if not l.startswith(("#", "!", "^"))]
g = pd.read_csv(StringIO("".join(gpl)), sep="\t", low_memory=False)
print(g.columns.tolist())
sym = g.set_index("ID").GENE_SYMBOL
mat = mat[mat.index.isin(sym.dropna().index)]
mat["sym"] = sym.reindex(mat.index).values
E = mat.groupby("sym").mean()
meta = pd.DataFrame({"gsm": gsm, "title": titles})
meta["patient"] = meta.title.str.extract(r"Patient ID_(\d+)")[0]
meta["tp"] = meta.title.str.extract(r"_(BL|C1D1|C1D15|Surg)$")[0]
print(meta.tp.value_counts())
rob = pd.read_csv(f"{OUT}/robust_DEGs.csv"); core = json.load(open(f"{OUT}/ml_selection.json"))["core"]
coef = pd.read_csv(f"{OUT}/risk_coefficients.csv", index_col=0).iloc[:, 0]
Z = E.sub(E.mean(1), axis=0).div(E.std(1), axis=0)
up = [x for x in rob.gene[rob.direction == "Up"] if x in Z.index]; dn = [x for x in rob.gene[rob.direction == "Down"] if x in Z.index]
meta = meta.set_index("gsm")
meta["resist_score"] = (Z.loc[up].mean() - Z.loc[dn].mean()).reindex(meta.index)
hub = [x for x in ["STAT1", "ISG15", "IFIT1", "IFIT2", "IFIT3", "IRF1", "IRF7", "CMPK2", "DDX60", "IFI44", "RSAD2", "OASL", "IFI6", "HELZ2"] if x in Z.index]
meta["ifn_hub_score"] = Z.loc[hub].mean().reindex(meta.index)
cg = [x for x in coef.index if x in Z.index]
meta["risk4"] = (Z.loc[cg].T @ coef[cg]).reindex(meta.index)
E2F = [l.split("\t") for l in open(f"{R}/MSigDB_Hallmark_2020.gmt") if l.startswith("E2F Targets")][0][2:]
E2F = [x.strip() for x in E2F if x.strip() in Z.index]
meta["e2f_score"] = Z.loc[E2F].mean().reindex(meta.index)
meta = meta.dropna(subset=["tp"]); meta.to_csv(f"{OUT}/neopalana_scores.csv")
out = {}
for a, b in [("C1D1", "C1D15"), ("BL", "C1D1"), ("C1D15", "Surg")]:
    w = meta[meta.tp.isin([a, b])].pivot_table(index="patient", columns="tp", values=["resist_score", "ifn_hub_score", "risk4", "e2f_score"])
    for sc in ["resist_score", "ifn_hub_score", "risk4", "e2f_score"]:
        x = w[sc].dropna()
        if len(x) < 5: continue
        st = stats.wilcoxon(x[a], x[b])
        out[f"{sc}:{a}->{b}"] = dict(n=len(x), median_change=float((x[b] - x[a]).median()), p=float(st.pvalue))
for k, v in out.items(): print(k, v)
json.dump(out, open(f"{OUT}/neopalana_results.json", "w"), indent=1)
print(len(up), len(dn), len(hub), len(E2F))
