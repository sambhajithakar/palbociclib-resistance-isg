"""Efficiency-corrected ddCq analysis for the interferon-signature qPCR validation.

Input
-----
plate_export.csv with columns:
    sample      biological sample name, e.g. "MCF7-PR_rep1"
    cell_model  "MCF7" or "T47D"
    condition   "parental" or "resistant"
    replicate   biological replicate index (1,2,3)
    gene        gene symbol
    well_type   "sample" | "NTC" | "noRT" | "std"
    dilution    for well_type == "std", the relative template amount (e.g. 1, 0.1, ...)
    Cq          quantification cycle from the instrument (blank/NaN if undetermined)

Everything else is derived. Run:  python ddct_analysis.py plate_export.csv out/

Method
------
1. QC gates: technical-replicate SD, NTC and -RT separation, Cq ceiling.
2. Amplification efficiency per assay from the standard-curve slope,
   E = 10^(-1/slope); accepted range 90-110 %.
3. Reference-gene stability by geNorm M (Vandesompele 2002) and NormFinder
   (Andersen 2004); the two most stable genes are used.
4. Efficiency-corrected relative quantity per Pfaffl (2001):
       Q(gene, sample) = E_gene ^ (Cq_calibrator_mean - Cq_sample)
   normalised by the geometric mean of the reference-gene quantities.
5. log2 fold-change resistant vs parental, one-sample t-test per cell model,
   Benjamini-Hochberg across target genes.

No data are supplied with this script; it is the analysis half of the protocol.
"""
import sys, os, math, itertools
import numpy as np, pandas as pd
from scipy import stats

TARGETS = ["PARP12", "CASP4", "SH3TC1", "LAMB1", "STAT1", "ISG15", "LGALS3BP", "MARCKS"]
REFS = ["ACTB", "GAPDH", "RPLP0", "TBP"]
CQ_MAX = 35.0
SD_MAX = 0.5


# ---------------------------------------------------------------- QC
def qc(df):
    log = []
    d = df.copy()
    d["Cq"] = pd.to_numeric(d["Cq"], errors="coerce")
    d.loc[d.Cq > CQ_MAX, "Cq"] = np.nan

    ntc = d[d.well_type == "NTC"].groupby("gene").Cq.min()
    nort = d[d.well_type == "noRT"].groupby(["gene", "sample"]).Cq.min()
    samp = d[d.well_type == "sample"]

    tech = samp.groupby(["gene", "sample"]).Cq.agg(["mean", "std", "count"])
    bad = tech[(tech["std"] > SD_MAX)]
    for (g, s), r in bad.iterrows():
        log.append(f"TECHNICAL SD  {g:9s} {s:16s} SD={r['std']:.2f} (>{SD_MAX}) - re-run")

    for g, c in ntc.items():
        if np.isfinite(c):
            worst = tech.xs(g, level=0)["mean"].max() if g in tech.index.get_level_values(0) else np.nan
            if np.isfinite(worst) and c < worst + 5:
                log.append(f"NTC           {g:9s} Cq={c:.1f} within 5 cycles of least abundant sample ({worst:.1f})")

    for (g, s), c in nort.items():
        if not np.isfinite(c):
            continue
        key = (g, s)
        if key in tech.index:
            plus = tech.loc[key, "mean"]
            if np.isfinite(plus) and c < plus + 5:
                log.append(f"-RT           {g:9s} {s:16s} Cq={c:.1f} vs +RT {plus:.1f} - genomic carry-over")
    return tech.reset_index(), log


# ------------------------------------------------------- efficiency
def efficiencies(df):
    out = []
    std = df[df.well_type == "std"].copy()
    std["Cq"] = pd.to_numeric(std["Cq"], errors="coerce")
    std["logq"] = np.log10(pd.to_numeric(std["dilution"], errors="coerce"))
    for g, s in std.groupby("gene"):
        s = s.dropna(subset=["Cq", "logq"])
        if len(s) < 6:
            out.append(dict(gene=g, slope=np.nan, E=np.nan, R2=np.nan, n=len(s), pass_=False))
            continue
        lr = stats.linregress(s.logq, s.Cq)
        E = 10 ** (-1 / lr.slope)
        eff_pct = (E - 1) * 100
        out.append(dict(gene=g, slope=lr.slope, E=E, efficiency_pct=eff_pct,
                        R2=lr.rvalue ** 2, n=len(s),
                        pass_=bool(90 <= eff_pct <= 110 and lr.rvalue ** 2 >= 0.99)))
    return pd.DataFrame(out)


# -------------------------------------------- reference-gene stability
def genorm_M(Q):
    """Q: samples x genes relative quantities. Returns M per gene (Vandesompele 2002)."""
    L = np.log2(Q)
    M = {}
    for g in Q.columns:
        vs = [np.std(L[g] - L[h], ddof=1) for h in Q.columns if h != g]
        M[g] = float(np.mean(vs))
    return pd.Series(M).sort_values()


def normfinder(Q, groups):
    """Simplified NormFinder stability (Andersen 2004): combines intra- and
    inter-group variation of log-transformed quantities."""
    L = np.log2(Q)
    res = {}
    for g in Q.columns:
        x = L[g]
        intra = np.mean([np.var(x[groups == k], ddof=1) for k in np.unique(groups)])
        means = np.array([np.mean(x[groups == k]) for k in np.unique(groups)])
        inter = np.var(means, ddof=1) if len(means) > 1 else 0.0
        res[g] = float(math.sqrt(intra) + abs(math.sqrt(max(inter, 0))))
    return pd.Series(res).sort_values()


# --------------------------------------------------------------- main
def main(path, outdir):
    os.makedirs(outdir, exist_ok=True)
    df = pd.read_csv(path)
    tech, log = qc(df)
    eff = efficiencies(df)
    eff.to_csv(f"{outdir}/assay_efficiency.csv", index=False)
    if not eff.pass_.all():
        for _, r in eff[~eff.pass_].iterrows():
            log.append(f"EFFICIENCY    {r.gene:9s} E={r.get('efficiency_pct', float('nan')):.1f}% "
                       f"R2={r.R2:.4f} - outside 90-110% / R2>=0.99")
    E = dict(zip(eff.gene, eff.E.fillna(2.0)))

    meta = df[df.well_type == "sample"].drop_duplicates("sample")[
        ["sample", "cell_model", "condition", "replicate"]].set_index("sample")
    cq = tech.pivot(index="sample", columns="gene", values="mean")
    cq = cq.reindex(meta.index)

    # relative quantity, calibrator = mean Cq of the parental samples of the same model
    Q = pd.DataFrame(index=cq.index, columns=cq.columns, dtype=float)
    for model, idx in meta.groupby("cell_model").groups.items():
        par = meta.loc[idx].query("condition == 'parental'").index
        for g in cq.columns:
            cal = cq.loc[par, g].mean()
            Q.loc[idx, g] = E.get(g, 2.0) ** (cal - cq.loc[idx, g])

    refQ = Q[[g for g in REFS if g in Q.columns]].dropna(axis=1, how="all")
    M = genorm_M(refQ)
    NF = normfinder(refQ, meta.loc[refQ.index, "condition"].values)
    rank = (M.rank() + NF.rank()).sort_values()
    chosen = list(rank.index[:2])
    pd.DataFrame({"geNorm_M": M, "NormFinder": NF, "combined_rank": rank}).to_csv(
        f"{outdir}/reference_gene_stability.csv")
    log.append(f"REFERENCE     normalising to geometric mean of {chosen[0]} and {chosen[1]}")

    norm = np.exp(np.log(Q[chosen]).mean(axis=1))
    NRQ = Q.div(norm, axis=0)
    NRQ.join(meta).to_csv(f"{outdir}/normalised_relative_quantities.csv")

    rows = []
    for model, idx in meta.groupby("cell_model").groups.items():
        m = meta.loc[idx]
        res_s = m.query("condition == 'resistant'").index
        par_s = m.query("condition == 'parental'").index
        for g in TARGETS:
            if g not in NRQ.columns:
                continue
            lr = np.log2(NRQ.loc[res_s, g].values)
            lp = np.log2(NRQ.loc[par_s, g].values)
            paired = min(len(lr), len(lp))
            if paired < 2:
                continue
            d = lr[:paired] - lp[:paired]
            t, p = stats.ttest_1samp(d, 0.0)
            rows.append(dict(cell_model=model, gene=g, n=paired,
                             log2FC=float(np.mean(d)), sd=float(np.std(d, ddof=1)),
                             t=float(t), p=float(p)))
    S = pd.DataFrame(rows)
    if len(S):
        S["FDR"] = np.nan
        for model, sub in S.groupby("cell_model"):
            o = np.argsort(sub.p.values)
            q = sub.p.values[o] * len(o) / (np.arange(len(o)) + 1)
            q = np.minimum.accumulate(q[::-1])[::-1]
            S.loc[sub.index[o], "FDR"] = np.clip(q, 0, 1)
        S = S.sort_values(["cell_model", "p"])
    S.to_csv(f"{outdir}/ddct_results.csv", index=False)

    with open(f"{outdir}/qc_report.txt", "w") as f:
        f.write("\n".join(log) if log else "All QC gates passed.\n")
    print("\n".join(log) if log else "All QC gates passed.")
    print()
    print(S.to_string(index=False) if len(S) else "no comparisons could be formed")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "qpcr_out")
