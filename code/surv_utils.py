"""Survival analysis in numpy/scipy: Cox PH (Breslow), LASSO-Cox, KM, log-rank, C-index, Uno time-AUC, DCA."""
import numpy as np, pandas as pd
from scipy import stats


def _risksets(time):
    o = np.argsort(-time, kind="mergesort")
    return o


def cox_fit(X, time, event, l2=0.0, max_iter=100, tol=1e-9):
    """Cox PH with Breslow ties via Newton-Raphson. Returns beta, se, loglik."""
    X = np.asarray(X, float); n, p = X.shape
    o = np.argsort(-time, kind="mergesort"); X, t, e = X[o], time[o], event[o]
    beta = np.zeros(p)
    def ll_grad_hess(b):
        eta = X @ b; eta -= eta.max(); w = np.exp(eta)
        S0 = np.cumsum(w); S1 = np.cumsum(w[:, None] * X, 0)
        S2 = np.cumsum(w[:, None, None] * X[:, :, None] * X[:, None, :], 0)
        # ties: use last index of each tied time block (all at-risk incl. ties)
        _, last = np.unique(t[::-1], return_index=True); last = len(t) - 1 - last
        idx = np.empty(n, int); ut = {tt: i for tt, i in zip(np.unique(t[::-1]), last)}
        for i in range(n): idx[i] = ut[t[i]]
        s0, s1, s2 = S0[idx], S1[idx], S2[idx]
        ll = np.sum(e * (eta - np.log(s0))) - 0.5 * l2 * b @ b
        m = s1 / s0[:, None]
        g = np.sum(e[:, None] * (X - m), 0) - l2 * b
        H = -np.sum(e[:, None, None] * (s2 / s0[:, None, None] - m[:, :, None] * m[:, None, :]), 0) - l2 * np.eye(p)
        return ll, g, H
    ll_old = -np.inf
    for _ in range(max_iter):
        ll, g, H = ll_grad_hess(beta)
        step = np.linalg.solve(H, g)
        nb = beta - step
        ll_new = ll_grad_hess(nb)[0]
        k = 0
        while ll_new < ll and k < 20:
            step /= 2; nb = beta - step; ll_new = ll_grad_hess(nb)[0]; k += 1
        beta = nb
        if abs(ll_new - ll) < tol: break
    ll, g, H = ll_grad_hess(beta)
    cov = np.linalg.inv(-H)
    return beta, np.sqrt(np.diag(cov)), ll, cov


def cox_summary(df, time, event, covars):
    d = df.dropna(subset=covars + [time, event])
    X = pd.get_dummies(d[covars], drop_first=True, dtype=float)
    b, se, ll, cov = cox_fit(X.values, d[time].values, d[event].values)
    z = b / se; p = 2 * stats.norm.sf(np.abs(z))
    return pd.DataFrame({"term": X.columns, "HR": np.exp(b), "lower": np.exp(b - 1.96 * se),
                         "upper": np.exp(b + 1.96 * se), "p": p, "coef": b}), ll, len(d), X


def cox_loglik(X, time, event, beta):
    o = np.argsort(-time, kind="mergesort"); X, t, e = X[o], time[o], event[o]
    eta = X @ beta; w = np.exp(eta - eta.max()); S0 = np.cumsum(w)
    # Breslow: at-risk set for tied times includes all ties
    ut, inv = np.unique(-t, return_inverse=True)
    last = np.zeros(len(ut), int)
    for i, k in enumerate(inv): last[k] = i
    s0 = S0[last[inv]]
    return np.sum(e * (eta - eta.max() - np.log(s0)))


def lasso_cox(X, time, event, lambdas=None, folds=10, seed=2026, max_iter=3000):
    """L1-penalised Cox by proximal gradient; lambda by K-fold CV of partial-likelihood deviance (Verweij–van Houwelingen)."""
    X = np.asarray(X, float); n, p = X.shape
    def fit(Xf, tf, ef, lam, b0=None):
        b = np.zeros(p) if b0 is None else b0.copy(); step = 1.0 / max(1, ef.sum())
        o = np.argsort(-tf, kind="mergesort"); Xo, eo = Xf[o], ef[o]
        for it in range(max_iter):
            eta = Xo @ b; w = np.exp(eta - eta.max()); S0 = np.cumsum(w); S1 = np.cumsum(w[:, None] * Xo, 0)
            g = -np.sum(eo[:, None] * (Xo - S1 / S0[:, None]), 0) / n
            bn = b - 0.5 * g
            bn = np.sign(bn) * np.maximum(np.abs(bn) - 0.5 * lam, 0)
            if np.max(np.abs(bn - b)) < 1e-7: b = bn; break
            b = bn
        return b
    if lambdas is None:
        o = np.argsort(-time); Xo, eo = X[o], event[o]
        S1 = np.cumsum(Xo, 0) / np.arange(1, n + 1)[:, None]
        g0 = np.abs(np.sum(eo[:, None] * (Xo - S1), 0) / n).max()
        lambdas = g0 * np.logspace(0, -2.5, 40)
    rng = np.random.default_rng(seed); fid = rng.permutation(np.arange(n) % folds)
    cvdev = np.zeros(len(lambdas))
    for k in range(folds):
        tr = fid != k; b = None
        for j, lam in enumerate(lambdas):
            b = fit(X[tr], time[tr], event[tr], lam, b)
            cvdev[j] += -2 * (cox_loglik(X, time, event, b) - cox_loglik(X[tr], time[tr], event[tr], b))
    best = int(np.argmin(cvdev)); b = None; path = []
    for lam in lambdas:
        b = fit(X, time, event, lam, b); path.append(b.copy())
    return np.array(path), lambdas, cvdev, best


def km(time, event):
    ut = np.unique(time[event == 1]); S = 1.0; out = [(0, 1.0, len(time))]
    for t in ut:
        n = (time >= t).sum(); d = ((time == t) & (event == 1)).sum()
        S *= 1 - d / n; out.append((t, S, n))
    return pd.DataFrame(out, columns=["time", "surv", "n_risk"])


def km_at(time, event, t0):
    k = km(time, event); return k.surv[k.time <= t0].iloc[-1]


def logrank(time, event, group):
    ut = np.unique(time[event == 1]); O = E = V = 0.0
    for t in ut:
        r = time >= t; n = r.sum(); n1 = (r & (group == 1)).sum()
        d = ((time == t) & (event == 1)).sum(); d1 = ((time == t) & (event == 1) & (group == 1)).sum()
        O += d1; E += d * n1 / n
        if n > 1: V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    chi = (O - E) ** 2 / V
    return chi, stats.chi2.sf(chi, 1)


def cindex(time, event, risk):
    num = den = 0.0
    for i in np.where(event == 1)[0]:
        m = time > time[i]
        den += m.sum(); num += (risk[i] > risk[m]).sum() + 0.5 * (risk[i] == risk[m]).sum()
    return num / den


def cindex_boot(time, event, risk, B=500, seed=2026):
    r = np.random.default_rng(seed); n = len(time); v = []
    for _ in range(B):
        i = r.integers(0, n, n); v.append(cindex(time[i], event[i], risk[i]))
    return np.percentile(v, [2.5, 97.5])


def _censor_surv(time, event):
    """KM of censoring distribution G(t)."""
    k = km(time, 1 - event)
    return lambda t: np.array([k.surv[k.time <= x].iloc[-1] if x >= 0 else 1.0 for x in np.atleast_1d(t)])


def time_auc(time, event, marker, t0):
    """Uno-type IPCW cumulative/dynamic AUC at t0 (Blanche et al. 2013, marginal weighting)."""
    G = _censor_surv(time, event)
    cases = (time <= t0) & (event == 1); ctrl = time > t0
    if cases.sum() == 0 or ctrl.sum() == 0: return np.nan
    w = 1 / np.maximum(G(time[cases] - 1e-9), 1e-6)
    mc, mk = marker[cases], marker[ctrl]
    comp = (mc[:, None] > mk[None, :]).astype(float) + 0.5 * (mc[:, None] == mk[None, :])
    return float((w[:, None] * comp).sum() / (w.sum() * ctrl.sum()))


def time_auc_ci(time, event, marker, t0, B=300, seed=2026):
    r = np.random.default_rng(seed); n = len(time); v = []
    for _ in range(B):
        i = r.integers(0, n, n); a = time_auc(time[i], event[i], marker[i], t0)
        if np.isfinite(a): v.append(a)
    return np.percentile(v, [2.5, 97.5])


def breslow_baseline(X, time, event, beta):
    eta = X @ beta; ut = np.unique(time[event == 1]); H = []; c = 0
    for t in ut:
        d = ((time == t) & (event == 1)).sum(); c += d / np.exp(eta[time >= t]).sum(); H.append(c)
    return ut, np.array(H)


def predict_risk_at(X, beta, ut, H, t0):
    H0 = H[ut <= t0][-1] if (ut <= t0).any() else 0.0
    return 1 - np.exp(-H0 * np.exp(X @ beta))


def net_benefit(time, event, prob, t0, thresholds):
    """Survival DCA (Vickers 2008): NB = TP/n - FP/n * pt/(1-pt), with KM-based event rates."""
    n = len(time); out = []
    Sall = km_at(time, event, t0)
    for pt in thresholds:
        pos = prob >= pt; npos = pos.sum()
        if npos == 0: nb = 0.0
        else:
            Spos = km_at(time[pos], event[pos], t0)
            tp = (1 - Spos) * npos / n; fp = Spos * npos / n
            nb = tp - fp * pt / (1 - pt)
        treat_all = (1 - Sall) - Sall * pt / (1 - pt)
        out.append((pt, nb, treat_all))
    return pd.DataFrame(out, columns=["threshold", "net_benefit", "treat_all"])
