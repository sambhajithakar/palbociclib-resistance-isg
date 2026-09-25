import numpy as np, pandas as pd, json
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.inspection import permutation_importance
from scipy import stats
OUT = "/home/claude/an/out"; SEED = 2026
rng = np.random.default_rng(SEED)
M = pd.read_pickle(f"{OUT}/merged_centered.pkl"); meta = pd.read_csv(f"{OUT}/meta.csv").set_index("sample")
cand = [g for g in open(f"{OUT}/ML_candidates.txt").read().split() if g in M.index]
X = M.loc[cand].T; X = (X - X.mean()) / X.std()
meta = meta.loc[X.index]; y = (meta.group == "Resistant").astype(int).values
print("samples", X.shape, "resistant", y.sum())
cv = StratifiedKFold(10, shuffle=True, random_state=SEED)

# (1) LASSO logistic
las = LogisticRegressionCV(Cs=np.logspace(-3, 2, 60), penalty="l1", solver="liblinear", cv=cv,
                           scoring="neg_log_loss", random_state=SEED, max_iter=5000).fit(X, y)
coef = pd.Series(las.coef_[0], index=cand)
lasso = list(coef[coef != 0].index)
# coefficient path for figure
Cs = np.logspace(-3, 2, 60); path = []
for c in Cs:
    path.append(LogisticRegression(penalty="l1", solver="liblinear", C=c, max_iter=5000).fit(X, y).coef_[0])
pd.DataFrame(path, index=Cs, columns=cand).to_csv(f"{OUT}/lasso_path.csv")
pd.DataFrame({"C": las.Cs_, "cv_logloss": -las.scores_[1].mean(0)}).to_csv(f"{OUT}/lasso_cv.csv", index=False)

# (2) SVM-RFE (Guyon 2002) with CV to choose subset size (ranking re-done inside each fold)
def svm_rank(Xa, ya):
    rem = list(range(Xa.shape[1])); order = []
    while rem:
        w = SVC(kernel="linear", C=1).fit(Xa[:, rem], ya).coef_[0]
        worst = rem[int(np.argmin(w ** 2))]; order.insert(0, worst); rem.remove(worst)
    return order
Xa = X.values; maxk = min(30, len(cand)); acc = np.zeros((10, maxk))
for f, (tr, te) in enumerate(cv.split(Xa, y)):
    rk = svm_rank(Xa[tr], y[tr])
    for k in range(1, maxk + 1):
        m = SVC(kernel="linear", C=1).fit(Xa[tr][:, rk[:k]], y[tr])
        acc[f, k - 1] = (m.predict(Xa[te][:, rk[:k]]) == y[te]).mean()
cvacc = acc.mean(0); bestk = int(np.argmax(cvacc)) + 1
svm = [cand[i] for i in svm_rank(Xa, y)[:bestk]]
pd.DataFrame({"k": range(1, maxk + 1), "cv_acc": cvacc}).to_csv(f"{OUT}/svmrfe_cv.csv", index=False)

# (3) Random forest + Boruta (Kursa & Rudnicki 2010)
def boruta(Xd, yv, iters=100):
    hits = np.zeros(Xd.shape[1]); r = np.random.default_rng(SEED)
    for i in range(iters):
        sh = np.apply_along_axis(r.permutation, 0, Xd)
        rf = RandomForestClassifier(500, random_state=int(r.integers(1e9)), n_jobs=-1).fit(np.hstack([Xd, sh]), yv)
        imp = rf.feature_importances_; p = Xd.shape[1]
        hits += imp[:p] > imp[p:].max()
    pvals = stats.binom.sf(hits - 1, iters, 0.5)
    padj = np.minimum(pvals * Xd.shape[1], 1)
    return hits, padj
hits, padj = boruta(Xa, y)
rf = RandomForestClassifier(1000, random_state=SEED, n_jobs=-1, oob_score=True).fit(Xa, y)
rfimp = pd.DataFrame({"gene": cand, "boruta_hits": hits, "boruta_padj": padj, "gini": rf.feature_importances_}).sort_values("gini", ascending=False)
rfimp.to_csv(f"{OUT}/rf_boruta.csv", index=False)
rfsel = list(rfimp.gene[rfimp.boruta_padj < 0.05])

sets = {"LASSO": lasso, "SVM-RFE": svm, "RF-Boruta": rfsel}
votes = pd.Series([g for s in sets.values() for g in s]).value_counts()
core = list(votes[votes == 3].index)
rule = "3/3"
if len(core) < 3: core = list(votes[votes >= 2].index); rule = ">=2/3"
print({k: len(v) for k, v in sets.items()}, "core", rule, core)
json.dump({"sets": sets, "core": core, "rule": rule, "bestk": bestk, "rf_oob": rf.oob_score_}, open(f"{OUT}/ml_selection.json", "w"), indent=1)

# performance: per gene AUC, combined logistic 10-fold CV, leave-one-dataset-out
def auc_ci(yt, p, B=2000):
    a = roc_auc_score(yt, p); bs = []
    r = np.random.default_rng(SEED)
    for _ in range(B):
        i = r.integers(0, len(yt), len(yt))
        if len(set(yt[i])) == 2: bs.append(roc_auc_score(yt[i], p[i]))
    return a, *np.percentile(bs, [2.5, 97.5])
rows = []
for g in core: rows.append((g, *auc_ci(y, X[g].values)))
pcv = cross_val_predict(LogisticRegression(max_iter=5000, C=1.0), X[core], y, cv=cv, method="predict_proba")[:, 1]
rows.append(("Combined (10-fold CV)", *auc_ci(y, pcv)))
lodo = np.zeros(len(y))
for d in meta.dataset.unique():
    te = (meta.dataset == d).values
    m = LogisticRegression(max_iter=5000).fit(X[core][~te], y[~te]); lodo[te] = m.predict_proba(X[core][te])[:, 1]
    rows.append((f"Leave-out {d}", *auc_ci(y[te], lodo[te])))
auc = pd.DataFrame(rows, columns=["model", "AUC", "CI_low", "CI_high"]); auc.to_csv(f"{OUT}/ml_auc.csv", index=False)
print(auc.to_string(index=False))
pd.DataFrame({"sample": X.index, "y": y, "p_cv": pcv, "p_lodo": lodo, "dataset": meta.dataset.values}).to_csv(f"{OUT}/ml_predictions.csv", index=False)

# (4) Gradient boosting + permutation importance (non-linear confirmation)
gb = GradientBoostingClassifier(n_estimators=300, max_depth=2, learning_rate=0.05, subsample=0.8, random_state=SEED)
pgb = cross_val_predict(gb, X, y, cv=cv, method="predict_proba")[:, 1]
pgb_core = cross_val_predict(gb, X[core], y, cv=cv, method="predict_proba")[:, 1]
gb.fit(X, y)
pi = permutation_importance(gb, X, y, n_repeats=50, random_state=SEED, scoring="roc_auc")
gbi = pd.DataFrame({"gene": cand, "perm_importance": pi.importances_mean, "sd": pi.importances_std}).sort_values("perm_importance", ascending=False)
gbi.to_csv(f"{OUT}/gbm_importance.csv", index=False)
print("GBM CV AUC all", roc_auc_score(y, pgb), "core", roc_auc_score(y, pgb_core))
print("core in GBM top10:", [g for g in gbi.gene.head(10) if g in core])
pd.DataFrame({"sample": X.index, "y": y, "gbm_all": pgb, "gbm_core": pgb_core}).to_csv(f"{OUT}/gbm_predictions.csv", index=False)
