"""Redocking validation of the grid definition: free phosphotyrosine into two candidate boxes."""
import pickle, numpy as np, pandas as pd
from vina import Vina
D = "/home/claude/an/dock"
lig = pickle.load(open(f"{D}/ligands.pkl", "rb")); ref = np.load(f"{D}/ref_ptyr.npy")
C = {"whole residue": np.load(f"{D}/center.npy"), "phosphate": np.load(f"{D}/phosphate_center.npy")}
rows = []
for box, key in [(24, "whole residue"), (20, "whole residue"), (18, "phosphate"), (16, "phosphate")]:
    v = Vina(sf_name="vina", cpu=2, seed=2026, verbosity=0); v.set_receptor(f"{D}/receptor.pdbqt")
    v.set_ligand_from_string(lig["phosphotyrosine (crystal ligand)"])
    v.compute_vina_maps(center=C[key].tolist(), box_size=[box] * 3)
    v.dock(exhaustiveness=16, n_poses=9)
    e = float(v.energies(n_poses=1)[0][0])
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in v.poses(n_poses=1).split("\n")
                    if l.startswith(("ATOM", "HETATM")) and l[77:].strip() != "HD"])
    d = np.linalg.norm(xyz[:, None, :] - ref[None, :, :], axis=2)
    rows.append(dict(box=box, centre=key, score=e, centroid_dist=float(np.linalg.norm(xyz.mean(0) - ref.mean(0))),
                     mean_nn_dist=float(d.min(1).mean()), max_nn_dist=float(d.min(1).max())))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f"{D}/redock_validation.csv", index=False)
