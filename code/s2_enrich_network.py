import numpy as np, pandas as pd, networkx as nx, math
from stats_utils import ora, gsea_prerank
R = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"; OUT = "/home/claude/an/out"

def gmt(f):
    d = {}
    for l in open(f"{R}/{f}"):
        p = l.rstrip("\n").split("\t")
        d[p[0]] = [x.split(",")[0] for x in p[2:] if x]
    return d
H = gmt("MSigDB_Hallmark_2020.gmt"); GO = gmt("GO_Biological_Process_2023.gmt"); KG = gmt("KEGG_2021_Human.gmt")
ra = pd.read_csv(f"{OUT}/RRA_all.csv"); rob = pd.read_csv(f"{OUT}/robust_DEGs.csv")
universe = sorted(set(ra.gene))
res = []
for d in ["Up", "Down"]:
    g = rob.gene[rob.direction == d]
    for lib, name in [(GO, "GO_BP"), (KG, "KEGG"), (H, "Hallmark")]:
        o = ora(g, lib, universe); o["direction"] = d; o["library"] = name; res.append(o)
E = pd.concat(res); E.to_csv(f"{OUT}/ORA_results.csv", index=False)
for (d, l), s in E.groupby(["direction", "library"]):
    print(f"\n== {d} {l}: sig FDR<0.05 = {(s.FDR<0.05).sum()}")
    print(s.head(8)[["term", "overlap", "set_size", "p", "FDR"]].to_string(index=False))

rank = ra.drop_duplicates("gene").set_index("gene").meanLFC.dropna()
gs = gsea_prerank(rank, H, nperm=2000)
gs.to_csv(f"{OUT}/GSEA_Hallmark.csv", index=False)
print("\n== GSEA Hallmark (FDR<0.05):"); print(gs[gs.FDR < 0.05][["term", "size", "NES", "p", "FDR"]].to_string(index=False))

# ---- network ----
net = pd.read_csv(f"{R}/STRING_network_robust.tsv", sep="\t")
G = nx.Graph()
for a, b, s in net[["preferredName_A", "preferredName_B", "score"]].itertuples(index=False):
    if s >= 0.4: G.add_edge(a, b, weight=s)
lcc = G.subgraph(max(nx.connected_components(G), key=len)).copy()
mcc = {v: 0.0 for v in lcc}
for c in nx.find_cliques(lcc):
    if len(c) >= 2:
        for v in c: mcc[v] += math.factorial(len(c) - 1)
C = pd.DataFrame({"Degree": dict(lcc.degree()), "MCC": mcc,
                  "Betweenness": nx.betweenness_centrality(lcc, normalized=True),
                  "Closeness": nx.closeness_centrality(lcc)})
tops = {m: set(C[m].sort_values(ascending=False).head(20).index) for m in C.columns}
C["n_top"] = [sum(g in t for t in tops.values()) for g in C.index]
C = C.sort_values(["n_top", "MCC"], ascending=False)
C["direction"] = rob.set_index("gene").direction.reindex(C.index).values
C["meanLFC"] = rob.set_index("gene").meanLFC.reindex(C.index).values
C.to_csv(f"{OUT}/network_centrality.csv")
hubs = list(C.index[C.n_top >= 3])
print(f"\nNetwork: all nodes {G.number_of_nodes()} edges {G.number_of_edges()}; LCC {lcc.number_of_nodes()} nodes {lcc.number_of_edges()} edges")
print("hubs:", len(hubs)); print(C.head(25).to_string())
cand = list(dict.fromkeys(hubs + list(rob.sort_values("score").gene.head(30))))
open(f"{OUT}/ML_candidates.txt", "w").write("\n".join(cand)); print("ML candidates", len(cand))
nx.write_gml(lcc, f"{OUT}/ppi_lcc.gml")
