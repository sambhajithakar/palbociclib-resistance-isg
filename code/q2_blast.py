"""In-silico PCR specificity check from NCBI BLAST (refseq_rna) hits.

A primer pair is judged specific when the only subject transcripts on which the
forward primer anneals to the plus strand and the reverse primer to the minus
strand, within a plausible product size, belong to the intended gene.
"""
import re, json, pandas as pd
R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"
OUT = "/home/claude/an/out"

prim = {}
name = None
for l in open("/home/claude/an/primers.fa"):
    l = l.strip()
    if l.startswith(">"): name = l[1:]
    elif l: prim[name] = l

rows = []
for l in open(f"{R}/blast_primers_tabular.txt"):
    if l.startswith("#") or "\t" not in l: continue
    p = l.rstrip("\n").split("\t")
    if len(p) < 12: continue
    rows.append(p)
H = pd.DataFrame(rows, columns=["q","s","pid","alen","mm","gap","qs","qe","ss","se","ev","bit"])
for c in ["pid","alen","mm","gap","qs","qe","ss","se","bit"]: H[c] = H[c].astype(float)
H["plen"] = H.q.map(lambda x: len(prim[x]))
print("total HSPs:", len(H), " queries:", H.q.nunique())

# an annealing-competent hit: >=95% identity over >= plen-1 nt, 3' end included
full = H[(H.pid >= 95) & (H.alen >= H.plen - 1) & (H.qe >= H.plen - 1)].copy()
full["strand"] = ["plus" if r.se > r.ss else "minus" for r in full.itertuples()]
print("annealing-competent hits:", len(full))

genes = sorted({q.rsplit("_", 1)[0] for q in prim})
pairs = []
for g in genes:
    F = full[full.q == g + "_F"]; Rv = full[full.q == g + "_R"]
    f_plus = F[F.strand == "plus"]; r_minus = Rv[Rv.strand == "minus"]
    common = set(f_plus.s) & set(r_minus.s)
    amps = []
    for acc in common:
        a = f_plus[f_plus.s == acc].ss.min(); b = r_minus[r_minus.s == acc].ss.max()
        size = b - a + 1
        if 50 <= size <= 2000: amps.append((acc, int(size)))
    pairs.append(dict(gene=g, n_F_sites=len(f_plus), n_R_sites=len(r_minus),
                      n_amplicons=len(amps),
                      accessions=";".join(sorted(a for a, _ in amps))))
P = pd.DataFrame(pairs)
P.to_csv(f"{OUT}/qpcr_insilico_pcr_raw.csv", index=False)
accs = sorted({a for row in P.accessions for a in row.split(";") if a} |
              {a for row in P.accessions.tolist() for a in row.split(";") if a})
accs = sorted({a for s in P.accessions if s for a in s.split(";")})
open(f"{OUT}/blast_accessions.txt", "w").write("\n".join(accs))
print(P.to_string(index=False))
print("\nunique amplifiable accessions to annotate:", len(accs))
