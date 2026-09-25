"""Stage 2: re-dock the best-scoring compounds and the controls at high exhaustiveness,
in triplicate, and enumerate protein-ligand contacts of the top pose."""
import os, pickle, numpy as np, pandas as pd
from vina import Vina
D = "/home/claude/an/dock"; os.makedirs(f"{D}/poses", exist_ok=True)
TOP, EXH, SEEDS, BOX = 20, 32, (1, 2, 3), 16
lig = pickle.load(open(f"{D}/ligands.pkl", "rb")); P = np.load(f"{D}/phosphate_center.npy").tolist()
props = pd.read_csv(f"{D}/library_properties.csv").set_index("name")
s1 = pd.read_csv(f"{D}/stage1_scores.csv").dropna(subset=["vina_kcal_mol"])
controls = ["phosphotyrosine (crystal ligand)", "palbociclib"]
hits = [c for c in s1.sort_values("vina_kcal_mol").compound if c not in controls][:TOP]
todo = controls + hits
print(f"stage 2 on {len(todo)} compounds", flush=True)

rec = [(l[17:20].strip() + l[22:26].strip(), l[12:16].strip(),
        np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l[77:].strip())
       for l in open(f"{D}/receptor.pdbqt") if l.startswith("ATOM")]
rxyz = np.array([r[2] for r in rec])

def dock(pdbqt, seed):
    v = Vina(sf_name="vina", cpu=2, seed=seed, verbosity=0)
    v.set_receptor(f"{D}/receptor.pdbqt"); v.set_ligand_from_string(pdbqt)
    v.compute_vina_maps(center=P, box_size=[BOX] * 3)
    v.dock(exhaustiveness=EXH, n_poses=5)
    return float(v.energies(n_poses=1)[0][0]), v.poses(n_poses=1)

import json
PART = f"{D}/stage2_partial.jsonl"
done = set()
if os.path.exists(PART):
    done = {json.loads(l)["compound"] for l in open(PART) if l.strip()}
print(f"{len(done)} already done", flush=True)
rows = [json.loads(l) for l in open(PART)] if os.path.exists(PART) else []
for n in [x for x in todo if x not in done]:
    es, pose = [], None
    for s in SEEDS:
        e, p = dock(lig[n], s); es.append(e)
        if s == SEEDS[0]: pose = p
    open(f"{D}/poses/{n.replace('/', '_').replace(' ', '_')[:50]}.pdbqt", "w").write(pose)
    at = [(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l[77:].strip())
          for l in pose.split("\n") if l.startswith(("ATOM", "HETATM"))]
    lxyz = np.array([a[0] for a in at]); ltyp = [a[1] for a in at]
    d = np.linalg.norm(lxyz[:, None, :] - rxyz[None, :, :], axis=2)
    close = {}
    for i, j in zip(*np.where(d < 4.0)):
        res = rec[j][0]; close[res] = min(close.get(res, 9), d[i, j])
    hb = set()
    POLAR = {"OA", "NA", "N", "HD", "SA"}
    for i, j in zip(*np.where(d < 3.5)):
        if ltyp[i] in POLAR and rec[j][3] in POLAR: hb.add(rec[j][0])
    p_ = props.loc[n] if n in props.index else None
    rows.append(dict(compound=n, vina_mean=np.mean(es), vina_sd=np.std(es, ddof=1),
                     vina_best=min(es), n_contacts=len(close),
                     contacts=";".join(sorted(close, key=lambda r: close[r])),
                     hbond_residues=";".join(sorted(hb)),
                     MW=None if p_ is None else p_.MW, LogP=None if p_ is None else p_.LogP,
                     QED=None if p_ is None else p_.QED, TPSA=None if p_ is None else p_.TPSA,
                     RotB=None if p_ is None else p_.RotB, HBD=None if p_ is None else p_.HBD,
                     HBA=None if p_ is None else p_.HBA,
                     Lipinski_viol=None if p_ is None else p_.Lipinski_viol,
                     smiles=None if p_ is None else p_.smiles, stage1=float(s1.set_index("compound").vina_kcal_mol.get(n, np.nan))))
    with open(PART, "a") as fh: fh.write(json.dumps(rows[-1], default=float) + "\n")
    print(f"{rows[-1]['vina_mean']:7.2f} ± {rows[-1]['vina_sd']:.2f}  {n[:45]}", flush=True)
R = pd.DataFrame(rows).sort_values("vina_mean")
ref = R.loc[R.compound == "phosphotyrosine (crystal ligand)", "vina_mean"]
if len(ref): R["better_than_pTyr"] = R.vina_mean < ref.values[0]
R.to_csv(f"{D}/stage2_top.csv", index=False)
print(R[["compound", "vina_mean", "vina_sd", "n_contacts", "hbond_residues"]].to_string(index=False))
