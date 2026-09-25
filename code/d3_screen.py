"""Stage 1 virtual screen of approved drugs against the validated STAT1 SH2 pTyr subsite.
Writes results incrementally and is resumable."""
import os, pickle, time, sys, numpy as np, pandas as pd
from vina import Vina
D = "/home/claude/an/dock"
BOX, EXH, SEED = 16, 4, 2026
lig = pickle.load(open(f"{D}/ligands.pkl", "rb"))
P = np.load(f"{D}/phosphate_center.npy").tolist()
props = pd.read_csv(f"{D}/library_properties.csv")

# pocket-matched subset: the native ligand is phosphotyrosine (MW 261, 16 heavy atoms)
keep = props[(props.MW <= 350) & (props.RotB <= 6) & (props.heavy <= 26)].name.tolist()
controls = ["phosphotyrosine (crystal ligand)", "palbociclib"]
todo = controls + [n for n in keep if n not in controls]
out = f"{D}/stage1_scores.csv"
done = set(pd.read_csv(out).compound) if os.path.exists(out) else set()
if not os.path.exists(out):
    open(out, "w").write("compound,vina_kcal_mol,seconds\n")
todo = [n for n in todo if n not in done and n in lig]
print(f"library {len(keep)} + controls; {len(done)} done, {len(todo)} to dock", flush=True)

t0 = time.time()
for i, n in enumerate(todo, 1):
    try:
        t = time.time()
        v = Vina(sf_name="vina", cpu=2, seed=SEED, verbosity=0)
        v.set_receptor(f"{D}/receptor.pdbqt")
        v.set_ligand_from_string(lig[n])
        v.compute_vina_maps(center=P, box_size=[BOX] * 3)
        v.dock(exhaustiveness=EXH, n_poses=3)
        e = float(v.energies(n_poses=1)[0][0]); dt = time.time() - t
    except Exception as ex:
        e, dt = float("nan"), 0.0
        print(f"  skip {n[:40]}: {str(ex)[:70]}", flush=True)
    with open(out, "a") as f:
        f.write('"%s",%.3f,%.1f\n' % (n.replace('"', "'"), e, dt))
    if i % 25 == 0:
        el = time.time() - t0
        print(f"{i}/{len(todo)}  {el/60:.0f} min elapsed, ~{el/i*(len(todo)-i)/3600:.1f} h left", flush=True)
print("stage 1 complete in %.2f h" % ((time.time() - t0) / 3600), flush=True)
