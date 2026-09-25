"""Prepare the STAT1 SH2 receptor (PDB 1BF5, chain A) and define the pTyr-binding grid box."""
import os, numpy as np, gemmi
from pdbfixer import PDBFixer
from openmm.app import PDBFile

RAW = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"
D = "/home/claude/an/dock"; os.makedirs(D, exist_ok=True)
POCKET = [584, 602, 603, 604, 605, 606, 607, 613, 630, 631, 632, 633]

# ---- box centre: pTyr701 of the partner protomer in the REMARK 350 biological dimer,
#      which occupies the SH2 pocket of chain A (Chen et al., Cell 1998).
st = gemmi.read_structure(f"{RAW}/1BF5.pdb"); st.setup_entities()
dimer = gemmi.Structure(); dimer.cell = st.cell; dimer.spacegroup_hm = st.spacegroup_hm
dimer.add_model(gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.AddNumber))
m2 = dimer[0]
pocket_xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in m2["A1"] if r.seqid.num in POCKET for a in r])
best = None
for c in m2:
    for r in c:
        if r.name != "PTR": continue
        p = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in r])
        d = np.linalg.norm(p.mean(0) - pocket_xyz.mean(0))
        if best is None or d < best[0]: best = (d, c.name, p)
dist, chain_id, ref_pos = best
center = ref_pos.mean(0).round(3)
print(f"partner pTyr701 ({chain_id}) centroid {center}, {dist:.2f} A from chain A SH2 pocket centroid")

# ---- receptor: chain A only, protonated at pH 7.4
fixer = PDBFixer(filename=f"{RAW}/1BF5.pdb")
fixer.removeChains([c.index for c in fixer.topology.chains() if c.id != "A"])
fixer.findMissingResidues(); fixer.missingResidues = {}      # do not model long absent loops
fixer.findNonstandardResidues(); fixer.replaceNonstandardResidues()   # PTR701 -> TYR
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open(f"{D}/receptor_H.pdb", "w"), keepIds=True)

# ---- AutoDock atom typing -> rigid receptor PDBQT (Vina scoring ignores partial charges)
AROM = {"PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"}, "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
        "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"}, "HIS": {"CG", "CD2", "CE1"}}
lines = [l for l in open(f"{D}/receptor_H.pdb") if l.startswith(("ATOM", "HETATM"))]
xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
elem = [l[76:78].strip() or l[12:16].strip()[0] for l in lines]
heavy_NOS = np.array([i for i, e in enumerate(elem) if e in ("N", "O", "S")])
out = []
for i, l in enumerate(lines):
    e, name, res = elem[i], l[12:16].strip(), l[17:20].strip()
    if e == "H":
        d = np.linalg.norm(xyz[heavy_NOS] - xyz[i], axis=1)
        if d.min() > 1.25: continue                       # non-polar H are merged away
        t = "HD"
    elif e == "C":
        t = "A" if name in AROM.get(res, ()) else "C"
    elif e == "N":
        if res == "HIS" and name in ("ND1", "NE2"):
            hd = [j for j, ej in enumerate(elem) if ej == "H" and np.linalg.norm(xyz[j] - xyz[i]) < 1.25]
            t = "N" if hd else "NA"
        else:
            t = "N"
    elif e == "O": t = "OA"
    elif e == "S": t = "SA"
    else: t = e
    out.append(f"{l[:54]}  1.00  0.00    {0.0:6.3f} {t:<2s}\n")
open(f"{D}/receptor.pdbqt", "w").writelines(["REMARK STAT1 SH2 domain, PDB 1BF5 chain A\n"] + out)
np.save(f"{D}/ref_ptyr.npy", ref_pos); np.save(f"{D}/center.npy", center)
print("receptor atoms written:", len(out), "of", len(lines))
print("pocket residues:", POCKET)
