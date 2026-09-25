"""Statistical helpers (Python re-implementations of limma eBayes, RRA, ORA, preranked GSEA)."""
import numpy as np
import pandas as pd
from scipy import special, stats


# ---------------------------------------------------------------- limma-style eBayes
def _trigamma_inverse(x):
    """Solve trigamma(y) = x (Smyth 2004, Newton iteration)."""
    x = np.asarray(x, float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(np.abs(dif / y) < 1e-8):
            break
    return y


def squeeze_var(s2, df):
    """Empirical Bayes prior (d0, s0^2) by method of moments on log variances (limma fitFDist)."""
    s2 = np.maximum(s2, 1e-12)
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = np.mean((e - emean) ** 2) * len(e) / (len(e) - 1) - special.polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s02 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf
        s02 = np.exp(emean)
    post = (d0 * s02 + df * s2) / (d0 + df) if np.isfinite(d0) else np.full_like(s2, s02)
    return d0, s02, post


def lm_ebayes(Y, design, coef):
    """Gene-wise linear model + empirical Bayes moderated t (limma lmFit + eBayes)."""
    X = np.asarray(design, float)
    n, p = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    B = Y @ X @ XtX_inv                       # genes x p
    resid = Y - B @ X.T
    df = n - p
    s2 = (resid ** 2).sum(1) / df
    d0, s02, s2post = squeeze_var(s2, df)
    se = np.sqrt(s2post * XtX_inv[coef, coef])
    t = B[:, coef] / se
    dft = df + (d0 if np.isfinite(d0) else 1e6)
    p = 2 * stats.t.sf(np.abs(t), dft)
    return B[:, coef], t, p, d0, s02


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def logcpm(counts, prior=1.0):
    lib = counts.sum(0)
    # TMM-free upper-quartile-robust: use median-of-ratios size factors (DESeq2 style)
    logc = np.log(counts.replace(0, np.nan))
    geo = logc.mean(1)
    ok = np.isfinite(geo)
    sf = np.exp((logc[ok].sub(geo[ok], axis=0)).median(0))
    norm = counts / sf
    eff_lib = lib.mean()
    return np.log2(norm / eff_lib * 1e6 + prior), sf


# ---------------------------------------------------------------- Robust rank aggregation
def rra(rank_lists, N):
    """Kolde et al. 2012 rho score with Bonferroni correction. rank_lists: list of ordered gene lists."""
    genes = sorted(set().union(*rank_lists))
    k = len(rank_lists)
    R = np.ones((len(genes), k))
    idx = {g: i for i, g in enumerate(genes)}
    for j, lst in enumerate(rank_lists):
        for r, g in enumerate(lst):
            R[idx[g], j] = (r + 1) / N
    R.sort(1)
    ks = np.arange(1, k + 1)
    pb = stats.beta.cdf(R, ks, k - ks + 1)
    rho = np.minimum(pb.min(1) * k, 1)
    return pd.DataFrame({"gene": genes, "score": rho}).sort_values("score").reset_index(drop=True)


# ---------------------------------------------------------------- ORA / GSEA
def ora(genes, gene_sets, universe, min_size=10, max_size=500):
    genes = set(genes) & set(universe)
    U = set(universe)
    N, n = len(U), len(genes)
    rows = []
    for name, gs in gene_sets.items():
        gs = set(gs) & U
        if not (min_size <= len(gs) <= max_size):
            continue
        k = len(gs & genes)
        if k == 0:
            continue
        p = stats.hypergeom.sf(k - 1, N, len(gs), n)
        rows.append((name, k, len(gs), n, p, ";".join(sorted(gs & genes))))
    df = pd.DataFrame(rows, columns=["term", "overlap", "set_size", "query_size", "p", "genes"])
    if len(df):
        df["FDR"] = bh(df.p.values)
        df["fold"] = (df.overlap / df.query_size) / (df.set_size / N)
    return df.sort_values("p").reset_index(drop=True)


def gsea_prerank(ranked: pd.Series, gene_sets, nperm=2000, min_size=15, max_size=500, seed=2026):
    """Preranked GSEA (Subramanian 2005), weighted KS (p=1), gene-set permutation."""
    rng = np.random.default_rng(seed)
    ranked = ranked.sort_values(ascending=False)
    g = ranked.index.values
    w = np.abs(ranked.values)
    N = len(g)
    pos = {x: i for i, x in enumerate(g)}

    def es(idx):
        hit = np.zeros(N, bool)
        hit[idx] = True
        ph = np.cumsum(np.where(hit, w, 0)) / w[hit].sum()
        pm = np.cumsum(~hit) / (N - hit.sum())
        d = ph - pm
        return d[np.argmax(np.abs(d))]

    rows = []
    for name, gs in gene_sets.items():
        idx = np.array([pos[x] for x in gs if x in pos])
        if not (min_size <= len(idx) <= max_size):
            continue
        obs = es(idx)
        null = np.array([es(rng.choice(N, len(idx), replace=False)) for _ in range(nperm)])
        same = null[np.sign(null) == np.sign(obs)]
        nes = obs / np.mean(np.abs(same)) if len(same) else np.nan
        p = (np.sum(np.abs(same) >= abs(obs)) + 1) / (len(same) + 1)
        rows.append((name, len(idx), obs, nes, p))
    df = pd.DataFrame(rows, columns=["term", "size", "ES", "NES", "p"])
    df["FDR"] = bh(df.p.values)
    return df.sort_values("p").reset_index(drop=True)
