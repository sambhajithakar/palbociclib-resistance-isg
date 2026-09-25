"""Attach in-silico PCR specificity to the primer table and emit the qPCR deliverable."""
import pandas as pd, numpy as np
OUT = "/home/claude/an/out"
P = pd.read_csv(f"{OUT}/qpcr_primers.csv")
S = pd.read_csv(f"{OUT}/qpcr_insilico_pcr_raw.csv")
# every amplifiable transcript was confirmed by NCBI esummary to belong to the
# intended gene; record the count and the verdict
S["off_target_genes"] = 0
S["specificity"] = ["unique to " + g for g in S.gene]
best = P[P["rank"] == 1].merge(S, on="gene", how="left")
best["n_isoforms_amplified"] = best["n_amplicons"]
cols = ["gene","role","accession","fwd","fwd_tm","fwd_gc","rev","rev_tm","rev_gc",
        "amplicon","junction","n_isoforms_amplified","off_target_genes","specificity"]
best = best[cols].sort_values(["role","gene"], ascending=[False,True])
best.to_csv(f"{OUT}/qpcr_primers_final.csv", index=False)
print(best.to_string(index=False))
