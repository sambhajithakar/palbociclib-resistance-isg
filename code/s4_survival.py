import numpy as np, pandas as pd, json
from scipy import stats
from surv_utils import *
R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"; OUT = "/home/claude/an/out"
core = json.load(open(f"{OUT}/ml_selection.json"))["core"]

# ---------------- TCGA ER+/HER2- ----------------
ex = pd.read_csv(f"{R}/TCGA_BRCA_HiSeqV2.gz", sep="\t", index_col=0)
cl = pd.read_csv(f"{R}/BRCA_clinicalMatrix.txt", sep="\t", low_memory=False).set_index("sampleID")
sv = pd.read_csv(f"{R}/BRCA_survival.txt", sep="\t").set_index("sample")
her2pos = (cl.lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive") | \
          (cl.lab_procedure_her2_neu_in_situ_hybrid_outcome_type == "Positive") | (cl.HER2_Final_Status_nature2012 == "Positive")
t = cl[(cl.breast_carcinoma_estrogen_receptor_status == "Positive") & ~her2pos & (cl.sample_type == "Primary Tumor")].copy()
t = t[t.index.isin(ex.columns) & t.index.isin(sv.index)]
t = t.join(sv[["OS", "OS.time", "PFI", "PFI.time"]])
t["PFI.time"] /= 365.25; t["OS.time"] /= 365.25
t = t[(t["PFI.time"] > 0) & t.PFI.notna()]
t["age"] = pd.to_numeric(t.age_at_initial_pathologic_diagnosis, errors="coerce")
st = t.pathologic_stage.fillna("")
t["stage"] = np.where(st.str.contains("IV|III"), "III-IV", np.where(st.str.contains("II"), "II", np.where(st.str.contains("Stage I"), "I", None)))
T = pd.DataFrame({"time": t["PFI.time"], "event": t.PFI.astype(int), "os_time": t["OS.time"], "os": t.OS,
                  "age": t.age, "stage": t.stage}, index=t.index)
print("TCGA ER+/HER2-:", len(T), "PFI events", T.event.sum(), "OS events", T.os.sum())

# ---------------- METABRIC ER+/HER2- ----------------
ml = pd.read_csv(f"{R}/METABRIC_mrna_long.tsv", sep="\t")
ann = pd.read_csv(f"{R}/Human.GRCh38.p13.annot.tsv.gz", sep="\t", usecols=["GeneID", "Symbol"]).set_index("GeneID").Symbol
ml["gene"] = ann.reindex(ml.entrez).values
mw = ml.pivot_table(index="gene", columns="sample", values="value", aggfunc="mean")
cp = pd.read_csv(f"{R}/METABRIC_clin_patient.tsv", sep="\t").pivot(index="patient", columns="attr", values="value")
cs = pd.read_csv(f"{R}/METABRIC_clin_sample.tsv", sep="\t")
pmap = cs.drop_duplicates("sample").set_index("sample").patient
cs = cs.pivot(index="sample", columns="attr", values="value")
m = cs.join(cp, on=pmap.reindex(cs.index).values) if False else cs.copy()
m["patient"] = pmap.reindex(m.index).values
m = m.join(cp, on="patient", rsuffix="_p")
m = m[(m.ER_STATUS == "Positive") & (m.HER2_STATUS == "Negative") & m.index.isin(mw.columns)]
M = pd.DataFrame({"time": pd.to_numeric(m.RFS_MONTHS, errors="coerce") / 12,
                  "event": m.RFS_STATUS.astype(str).str[0].map({"1": 1, "0": 0}),
                  "os_time": pd.to_numeric(m.OS_MONTHS, errors="coerce") / 12,
                  "os": m.OS_STATUS.astype(str).str[0].map({"1": 1, "0": 0}),
                  "age": pd.to_numeric(m.AGE_AT_DIAGNOSIS, errors="coerce"),
                  "stage": pd.to_numeric(m.TUMOR_STAGE, errors="coerce").map({0: "I", 1: "I", 2: "II", 3: "III-IV", 4: "III-IV"}),
                  "grade": pd.to_numeric(m.GRADE, errors="coerce"),
                  "hormone": m.HORMONE_THERAPY, "chemo": m.CHEMOTHERAPY}, index=m.index)
M = M[(M.time > 0) & M.event.notna()]; M["event"] = M.event.astype(int)
print("METABRIC ER+/HER2-:", len(M), "RFS events", M.event.sum(), "OS events", M.os.sum())

def z(expr, genes, samples):
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g, samples].T.astype(float)
    return (sub - sub.mean()) / sub.std()
Zt = z(ex, core, T.index); Zm = z(mw, core, M.index)
genes = [g for g in core if g in Zt.columns and g in Zm.columns]
print("core genes available in both cohorts:", len(genes), "missing:", set(core) - set(genes))
Zt, Zm = Zt[genes], Zm[genes]

# univariable Cox per gene (TCGA and METABRIC)
uni = []
for g in genes:
    for nm, Z, D in [("TCGA", Zt, T), ("METABRIC", Zm, M)]:
        b, se, ll, _ = cox_fit(Z[[g]].values, D.time.values, D.event.values)
        uni.append((g, nm, np.exp(b[0]), np.exp(b[0] - 1.96 * se[0]), np.exp(b[0] + 1.96 * se[0]), 2 * stats.norm.sf(abs(b[0] / se[0]))))
uni = pd.DataFrame(uni, columns=["gene", "cohort", "HR", "lower", "upper", "p"]); uni.to_csv(f"{OUT}/univariable_gene_cox.csv", index=False)
print(uni.pivot(index="gene", columns="cohort", values=["HR", "p"]).round(3))

# LASSO-Cox in TCGA
path, lams, cvdev, best = lasso_cox(Zt.values, T.time.values, T.event.values)
coef = pd.Series(path[best], index=genes)
sel = coef[coef != 0]
if len(sel) == 0:
    j = int(np.argmax([np.sum(p != 0) >= 3 for p in path])); coef = pd.Series(path[j], index=genes); sel = coef[coef != 0]
    print("lambda.min kept 0 genes; using first lambda with >=3 genes")
# refit unpenalised Cox on selected genes for final coefficients? keep LASSO coefficients (as in glmnet practice)
print("LASSO-Cox selected:", sel.round(4).to_dict())
pd.DataFrame(path, columns=genes, index=lams).to_csv(f"{OUT}/lassocox_path.csv")
pd.DataFrame({"lambda": lams, "cv_deviance": cvdev}).to_csv(f"{OUT}/lassocox_cv.csv", index=False)
sel.to_csv(f"{OUT}/risk_coefficients.csv")
T["risk"] = Zt[sel.index].values @ sel.values; M["risk"] = Zm[sel.index].values @ sel.values
for D in (T, M): D["high"] = (D.risk > D.risk.median()).astype(int)

res = {"n_tcga": len(T), "ev_tcga": int(T.event.sum()), "n_mb": len(M), "ev_mb": int(M.event.sum()), "coef": sel.to_dict()}
for nm, D, tt in [("TCGA", T, [3, 5, 8]), ("METABRIC", M, [3, 5, 10])]:
    for ep, tcol, ecol in [("primary", "time", "event"), ("OS", "os_time", "os")]:
        d = D.dropna(subset=[tcol, ecol]); tm, ev = d[tcol].values, d[ecol].values.astype(int)
        chi, p = logrank(tm, ev, d.high.values)
        b, se, _, _ = cox_fit(d[["high"]].values, tm, ev)
        bc, sec, _, _ = cox_fit(((d.risk - d.risk.mean()) / d.risk.std()).values[:, None], tm, ev)
        c = cindex(tm, ev, d.risk.values); cci = cindex_boot(tm, ev, d.risk.values, B=200)
        res[f"{nm}_{ep}"] = dict(n=len(d), events=int(ev.sum()), logrank_p=p, HR_group=np.exp(b[0]), HR_lo=np.exp(b[0]-1.96*se[0]),
                                 HR_hi=np.exp(b[0]+1.96*se[0]), HR_perSD=np.exp(bc[0]), HRsd_lo=np.exp(bc[0]-1.96*sec[0]),
                                 HRsd_hi=np.exp(bc[0]+1.96*sec[0]), p_perSD=2*stats.norm.sf(abs(bc[0]/sec[0])), C=c, C_lo=cci[0], C_hi=cci[1])
        if ep == "primary":
            res[f"{nm}_tAUC"] = {y: [time_auc(tm, ev, d.risk.values, y), *time_auc_ci(tm, ev, d.risk.values, y, B=200)] for y in tt}
        print(nm, ep, {k: round(v, 4) if isinstance(v, float) else v for k, v in res[f"{nm}_{ep}"].items()})
    print(nm, "tAUC", res[f"{nm}_tAUC"])

# multivariable Cox
mv = []
Tc = T.copy(); Tc["risk_sd"] = (Tc.risk - Tc.risk.mean()) / Tc.risk.std()
Mc = M.copy(); Mc["risk_sd"] = (Mc.risk - Mc.risk.mean()) / Mc.risk.std(); Mc["grade"] = Mc.grade.map({1: "1-2", 2: "1-2", 3: "3"})
for nm, D, cov in [("TCGA", Tc, ["risk_sd", "age", "stage"]), ("METABRIC", Mc, ["risk_sd", "age", "stage", "grade"])]:
    for v in cov:
        s, _, n, _ = cox_summary(D, "time", "event", [v]); s["model"] = "Univariable"; s["cohort"] = nm; mv.append(s)
    s, _, n, _ = cox_summary(D, "time", "event", cov); s["model"] = "Multivariable"; s["cohort"] = nm; s["n"] = n; mv.append(s)
mv = pd.concat(mv); mv.to_csv(f"{OUT}/cox_uni_multi.csv", index=False); print(mv.round(4).to_string(index=False))

# nomogram model, C-index comparison, LRT, CV-calibration and DCA in METABRIC
d = Mc.dropna(subset=["age", "stage", "grade"]).copy()
Xc = pd.get_dummies(d[["age", "stage", "grade"]], drop_first=True, dtype=float); Xf = Xc.copy(); Xf.insert(0, "risk_sd", d.risk_sd)
bc, _, llc, _ = cox_fit(Xc.values, d.time.values, d.event.values); bf, sef, llf, _ = cox_fit(Xf.values, d.time.values, d.event.values)
lrt = 2 * (llf - llc); p_lrt = stats.chi2.sf(lrt, 1)
Cc = cindex(d.time.values, d.event.values, Xc.values @ bc); Cf = cindex(d.time.values, d.event.values, Xf.values @ bf)
res["nomogram"] = dict(n=len(d), events=int(d.event.sum()), C_clin=Cc, C_full=Cf, LRT=lrt, p_LRT=p_lrt,
                       coef=dict(zip(Xf.columns, bf)), se=dict(zip(Xf.columns, sef)))
print("nomogram", res["nomogram"])
# 10-fold CV predicted 5-y risk
rng = np.random.default_rng(2026); fid = rng.permutation(np.arange(len(d)) % 10)
pc = np.zeros(len(d)); pf = np.zeros(len(d))
for k in range(10):
    tr, te = fid != k, fid == k
    for X_, arr in [(Xc.values, pc), (Xf.values, pf)]:
        b, _, _, _ = cox_fit(X_[tr], d.time.values[tr], d.event.values[tr])
        ut, H = breslow_baseline(X_[tr], d.time.values[tr], d.event.values[tr], b)
        arr[te] = predict_risk_at(X_[te], b, ut, H, 5)
d["p5_clin"], d["p5_full"] = pc, pf
q = pd.qcut(d.p5_full, 5, labels=False); cal = []
for g in range(5):
    s = d[q == g]; cal.append((s.p5_full.mean(), 1 - km_at(s.time.values, s.event.values, 5), len(s)))
cal = pd.DataFrame(cal, columns=["predicted", "observed", "n"]); cal.to_csv(f"{OUT}/calibration_5y.csv", index=False)
thr = np.round(np.arange(0.02, 0.41, 0.01), 2)
dca = net_benefit(d.time.values, d.event.values, d.p5_full.values, 5, thr).rename(columns={"net_benefit": "full"})
dca["clinical"] = net_benefit(d.time.values, d.event.values, d.p5_clin.values, 5, thr).net_benefit
dca.to_csv(f"{OUT}/dca_5y.csv", index=False)
print(cal); print(dca.iloc[::5])
T.to_csv(f"{OUT}/tcga_cohort.csv"); M.to_csv(f"{OUT}/metabric_cohort.csv"); d.to_csv(f"{OUT}/metabric_nomogram_data.csv")
json.dump(res, open(f"{OUT}/survival_results.json", "w"), indent=1, default=float)
