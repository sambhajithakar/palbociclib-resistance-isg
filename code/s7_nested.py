"""Fully nested leave-one-dataset-out validation: differential expression, rank aggregation and
all three feature selectors are re-run inside each training fold, so the held-out dataset is
never seen during feature selection."""
import numpy as np, pandas as pd, json
from scipy import stats
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from stats_utils import lm_ebayes, bh, rra
O = "/home/claude/an/out"; SEED = 2026
expr = pd.read_pickle(f"{O}/expr_all.pkl")            # per-dataset logCPM
meta = pd.read_csv(f"{O}/meta.csv").set_index("sample")
DS = list(expr.keys())
universe = sorted(set.intersection(*[set(e.index) for e in expr.values()]))

def de(d):
    m = meta[meta.dataset == d]
    y = m.loc[expr[d].columns.intersection(m.index)]
    Y = expr[d][y.index].reindex(universe).dropna()
    grp = (y.group == "Resistant").astype(int).values
    X = [np.ones(len(grp)), grp]
    if y.cell.nunique() > 1 and d == "GSE222367": X.append((y.cell == "T47D").astype(int).values)
    lfc, t, p, _, _ = lm_ebayes(Y.values, np.column_stack(X), 1)
    return pd.DataFrame({"gene": Y.index, "log2FC": lfc, "P": p}).sort_values("P")

DE = {d: de(d) for d in DS}

def boruta(Xd, yv, iters=60, seed=SEED):
    r = np.random.default_rng(seed); hits = np.zeros(Xd.shape[1])
    for _ in range(iters):
        sh = np.apply_along_axis(r.permutation, 0, Xd)
        rf = RandomForestClassifier(300, random_state=int(r.integers(1e9)), n_jobs=-1).fit(np.hstack([Xd, sh]), yv)
        imp = rf.feature_importances_; p = Xd.shape[1]
        hits += imp[:p] > imp[p:].max()
    return np.minimum(stats.binom.sf(hits - 1, iters, 0.5) * Xd.shape[1], 1)

def svm_rfe(Xa, ya, kmax=20):
    rem = list(range(Xa.shape[1])); order = []
    while rem:
        w = SVC(kernel="linear", C=1).fit(Xa[:, rem], ya).coef_[0]
        worst = rem[int(np.argmin(w ** 2))]; order.insert(0, worst); rem.remove(worst)
    cv = StratifiedKFold(5, shuffle=True, random_state=SEED); acc = []
    for k in range(1, min(kmax, Xa.shape[1]) + 1):
        a = []
        for tr, te in cv.split(Xa, ya):
            a.append((SVC(kernel="linear", C=1).fit(Xa[tr][:, order[:k]], ya[tr]).predict(Xa[te][:, order[:k]]) == ya[te]).mean())
        acc.append(np.mean(a))
    return order[:int(np.argmax(acc)) + 1]

def centred(d, genes):
    m = meta[meta.dataset == d]
    y = m.loc[expr[d].columns.intersection(m.index)]
    E = expr[d][y.index].reindex(genes)
    out = []
    for b, mm in (y.groupby("cell") if d != "GSE229235" else [("all", y)]):
        s = E[mm.index]; out.append(s.sub(s.mean(1), axis=0))
    return pd.concat(out, axis=1), (y.group == "Resistant").astype(int).values, y

rows = []
for held in DS:
    tr_ds = [d for d in DS if d != held]
    up = rra([list(DE[d][DE[d].log2FC > 0].gene) for d in tr_ds], len(universe))
    dn = rra([list(DE[d][DE[d].log2FC < 0].gene) for d in tr_ds], len(universe))
    lf = pd.DataFrame({d: DE[d].set_index("gene").log2FC.reindex(universe) for d in tr_ds})
    both = pd.concat([up.assign(direction="Up"), dn.assign(direction="Down")])
    both["meanLFC"] = lf.mean(1).reindex(both.gene).values
    both["ok"] = np.sign(lf).reindex(both.gene).apply(lambda s: abs(s.sum()) == len(s), axis=1).values
    rob = both[(both.score < 0.05) & (both.meanLFC.abs() > 1) & both.ok]
    cand = list(rob.sort_values("score").gene.head(40))

    Xtr, ytr, mtr = [], [], []
    for d in tr_ds:
        Xd, yd, md = centred(d, cand); Xtr.append(Xd); ytr.append(yd); mtr.append(md)
    Xtr = pd.concat(Xtr, axis=1).dropna(); ytr = np.concatenate(ytr)
    keep = Xtr.index.tolist()
    Z = Xtr.T; Z = (Z - Z.mean()) / Z.std()
    las = LogisticRegressionCV(Cs=np.logspace(-3, 2, 40), penalty="l1", solver="liblinear",
                               cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                               scoring="neg_log_loss", random_state=SEED, max_iter=5000).fit(Z, ytr)
    s_las = set(np.array(keep)[las.coef_[0] != 0])
    s_svm = set(np.array(keep)[svm_rfe(Z.values, ytr)])
    s_rf = set(np.array(keep)[boruta(Z.values, ytr) < 0.05])
    votes = pd.Series([g for s in (s_las, s_svm, s_rf) for g in s]).value_counts()
    fold_core = list(votes[votes >= 2].index) or list(votes.index[:5])

    Zc = Z[fold_core]
    clf = LogisticRegression(max_iter=5000).fit(Zc, ytr)
    Xte, yte, _ = centred(held, fold_core)
    Zt = Xte.T; Zt = (Zt - Zt.mean()) / Zt.std()
    p = clf.predict_proba(Zt[fold_core])[:, 1]
    a = roc_auc_score(yte, p)
    rows.append(dict(held_out=held, n_test=len(yte), n_robust=len(rob), n_core=len(fold_core),
                     core=";".join(sorted(fold_core)), AUC=a,
                     overlap_with_final=len(set(fold_core) & set(json.load(open(f"{O}/ml_selection.json"))["core"]))))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows); R.to_csv(f"{O}/nested_lodo.csv", index=False)
print(R[["held_out", "n_test", "n_robust", "n_core", "AUC", "overlap_with_final"]].round(3).to_string(index=False))
print("mean nested AUC %.3f" % R.AUC.mean())
