"""Design exon-junction-spanning qPCR primers for the core signature and reference genes.

Thermodynamics: SantaLucia (1998) unified nearest-neighbour parameters, salt-corrected
(Owczarzy 2004) for 50 mM monovalent, 3 mM Mg2+, 0.2 uM primer.
Constraints follow standard SYBR-Green qPCR practice; final specificity must still be
confirmed in NCBI Primer-BLAST before ordering.
"""
import re, math, json, itertools, numpy as np, pandas as pd
R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"
OUT = "/home/claude/an/out"

# ---- SantaLucia 1998 unified NN parameters: dH (kcal/mol), dS (cal/mol/K)
NN = {"AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3), "CA": (-8.5, -22.7),
      "GT": (-8.4, -22.4), "CT": (-7.8, -21.0), "GA": (-8.2, -22.2), "CG": (-10.6, -27.2),
      "GC": (-9.8, -24.4), "GG": (-8.0, -19.9)}
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def _nn(d):
    if d in NN: return NN[d]
    return NN[rc(d)]
def tm(seq, Ct=0.2e-6, Na_mM=50.0, Mg_mM=3.0, dNTP_mM=0.8):
    s = seq.upper()
    dH = 0.0; dS = 0.0
    for i in range(len(s) - 1):
        h, e = _nn(s[i:i+2]); dH += h; dS += e
    # initiation
    dH += 0.1 if s[0] in "GC" else 2.3
    dS += -2.8 if s[0] in "GC" else 4.1
    dH += 0.1 if s[-1] in "GC" else 2.3
    dS += -2.8 if s[-1] in "GC" else 4.1
    # equivalent monovalent concentration from Mg2+ (von Ahsen 2001), in mol/L
    na_eq = (Na_mM + 120 * math.sqrt(max(Mg_mM - dNTP_mM, 0.0))) / 1000.0
    t = (dH * 1000) / (dS + 1.987 * math.log(Ct / 4)) - 273.15
    gc = (s.count("G") + s.count("C")) / len(s)
    # Owczarzy 2004 salt correction
    inv = 1 / (t + 273.15) + ((4.29 * gc - 3.95) * 1e-5 * math.log(na_eq)
                              + 9.40e-6 * math.log(na_eq) ** 2)
    return 1 / inv - 273.15
def gc_frac(s): return (s.count("G") + s.count("C")) / len(s)
def max_run(s): return max(len(m.group(0)) for m in re.finditer(r"(.)\1*", s))
def self_comp(a, b, min_match=5):
    """Longest complementary stretch between a and reverse of b (crude dimer check)."""
    best = 0; b_rc = rc(b)
    for i in range(len(a)):
        for j in range(len(b_rc)):
            k = 0
            while i + k < len(a) and j + k < len(b_rc) and a[i+k] == b_rc[j+k]: k += 1
            best = max(best, k)
    return best
def three_prime_dimer(a, b, n=5):
    return self_comp(a[-n:], b[-n:], 3)

# ---- parse RefSeq records
recs = {}
for block in open(f"{R}/refseq_records.txt").read().split(">>>GENE ")[1:]:
    gene = block.split("\n")[0].strip()
    acc = re.search(r"^VERSION\s+(\S+)", block, re.M).group(1)
    exons = [tuple(int(x) for x in m) for m in
             re.findall(r"^\s{5}exon\s+<?(\d+)\.\.>?(\d+)", block, re.M)]
    cds = re.search(r"^\s{5}CDS\s+<?(\d+)\.\.>?(\d+)", block, re.M)
    seq = "".join(re.findall(r"^\s*\d+\s([acgtn ]+)$", block, re.M)).replace(" ", "").upper()
    recs[gene] = dict(acc=acc, exons=exons, seq=seq,
                      cds=(int(cds.group(1)), int(cds.group(2))) if cds else None)
    print(f"{gene:10s} {acc:16s} len={len(seq):6d} exons={len(exons)}")

def junctions(exons):
    return [e[1] for e in exons[:-1]]           # 1-based last base of each exon

def design(gene, amp=(80, 160), tm_range=(58.5, 62.0), plen=(19, 24), max_pairs=3, require_junction=True):
    r = recs[gene]; s = r["seq"]; J = junctions(r["exons"])
    if not s: return []
    lo, hi = (r["cds"] if r["cds"] else (1, len(s)))
    cands_f, cands_r = [], []
    for start in range(lo - 1, min(hi, len(s)) - plen[0]):
        for L in range(plen[0], plen[1] + 1):
            p = s[start:start+L]
            if len(p) < L or "N" in p: continue
            if not (0.40 <= gc_frac(p) <= 0.62): continue
            if max_run(p) > 3: continue
            if p[-1] not in "GC": continue
            if sum(c in "GC" for c in p[-5:]) > 3: continue
            t = tm(p)
            if not (tm_range[0] <= t <= tm_range[1]): continue
            if self_comp(p, p) > 6: continue
            cands_f.append((start + 1, start + L, p, t))
    for end in range(lo + plen[0], min(hi, len(s)) + 1):
        for L in range(plen[0], plen[1] + 1):
            st = end - L
            if st < 0: continue
            sub = s[st:end]
            if "N" in sub: continue
            p = rc(sub)
            if not (0.40 <= gc_frac(p) <= 0.62): continue
            if max_run(p) > 3: continue
            if p[-1] not in "GC": continue
            if sum(c in "GC" for c in p[-5:]) > 3: continue
            t = tm(p)
            if not (tm_range[0] <= t <= tm_range[1]): continue
            if self_comp(p, p) > 6: continue
            cands_r.append((st + 1, end, p, t))
    out = []
    for f in cands_f:
        for rv in cands_r:
            size = rv[1] - f[0] + 1
            if not (amp[0] <= size <= amp[1]): continue
            if abs(f[3] - rv[3]) > 1.5: continue
            if three_prime_dimer(f[2], rv[2]) >= 4: continue
            if self_comp(f[2], rv[2]) > 7: continue
            span_f = any(f[0] <= j < f[1] and (j - f[0] + 1) >= 4 and (f[1] - j) >= 4 for j in J)
            span_r = any(rv[0] <= j < rv[1] and (j - rv[0] + 1) >= 4 and (rv[1] - j) >= 4 for j in J)
            cross = any(f[1] <= j < rv[0] for j in J)
            if require_junction and not (span_f or span_r or cross): continue
            score = (2 if (span_f or span_r) else 0) + (1 if cross else 0) \
                    - abs(f[3] - rv[3]) - abs(gc_frac(f[2]) - 0.5) - abs(gc_frac(rv[2]) - 0.5)
            out.append(dict(gene=gene, accession=r["acc"], fwd=f[2], fwd_tm=round(f[3], 1),
                            fwd_gc=round(gc_frac(f[2]) * 100), fwd_pos=f[0],
                            rev=rv[2], rev_tm=round(rv[3], 1), rev_gc=round(gc_frac(rv[2]) * 100),
                            rev_pos=rv[1], amplicon=size,
                            junction=("primer spans junction" if (span_f or span_r)
                                      else "amplicon spans junction" if cross
                                      else "NO junction - use -RT control"),
                            score=round(score, 3)))
    out.sort(key=lambda d: -d["score"])
    # keep non-overlapping alternatives
    keep = []
    for o in out:
        if all(abs(o["fwd_pos"] - k["fwd_pos"]) > 30 for k in keep): keep.append(o)
        if len(keep) >= max_pairs: break
    return keep

TARGETS = ["PARP12", "CASP4", "SH3TC1", "LAMB1", "STAT1", "ISG15", "LGALS3BP", "MARCKS"]
REFS = ["ACTB", "GAPDH", "RPLP0", "TBP"]
rows = []
for g in TARGETS + REFS:
    d = design(g)
    if not d: d = design(g, amp=(70, 200), tm_range=(57.5, 63.0), plen=(18, 25))
    if not d: d = design(g, amp=(70, 220), tm_range=(57.0, 64.0), plen=(18, 26))
    if not d: d = design(g, amp=(70, 220), tm_range=(57.0, 64.0), plen=(18, 26), require_junction=False)
    for i, x in enumerate(d):
        x["role"] = "target" if g in TARGETS else "reference"
        x["rank"] = i + 1
        rows.append(x)
    print(f"{g:10s} {len(d)} pair(s)" + ("" if d else "   <-- none found"))
P = pd.DataFrame(rows)
P.to_csv(f"{OUT}/qpcr_primers.csv", index=False)
print("\nBest pair per gene:")
print(P[P["rank"] == 1][["gene", "accession", "fwd", "fwd_tm", "rev", "rev_tm", "amplicon", "junction"]].to_string(index=False))
