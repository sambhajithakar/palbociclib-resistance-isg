"""Build the docking library: approved small-molecule drugs (ChEMBL) + phosphotyrosine control."""
import pickle, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from meeko import MoleculePreparation, PDBQTWriterLegacy
RDLogger.DisableLog("rdApp.*")
RAW = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw"; D = "/home/claude/an/dock"

df = pd.read_csv(f"{RAW}/chembl_approved.tsv", sep="\t")
df = df[df.withdrawn.fillna(0).astype(float) == 0]
print("approved small molecules:", len(df))
chooser = rdMolStandardize.LargestFragmentChooser()
uncharger = rdMolStandardize.Uncharger()
prep = MoleculePreparation()

CONTROLS = {"phosphotyrosine (crystal ligand)": "N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O",
            "palbociclib": "CC(=O)c1c(C)c2cnc(Nc3ccc(N4CCNCC4)cn3)nc2n(C2CCCC2)c1=O"}

rows, lig = [], {}
def add(name, smi, src):
    m = Chem.MolFromSmiles(smi)
    if m is None: return
    m = uncharger.uncharge(chooser.choose(m))
    na = m.GetNumHeavyAtoms()
    if not (10 <= na <= 50): return
    mw = Descriptors.MolWt(m)
    if mw > 700: return
    mh = Chem.AddHs(m)
    p = AllChem.ETKDGv3(); p.randomSeed = 2026
    if AllChem.EmbedMolecule(mh, p) != 0: return
    try: AllChem.MMFFOptimizeMolecule(mh, maxIters=600)
    except Exception: pass
    try:
        setups = prep.prepare(mh)
        s, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    except Exception: return
    if not ok: return
    lig[name] = s
    rows.append(dict(name=name, source=src, smiles=Chem.MolToSmiles(m), MW=mw,
                     LogP=Descriptors.MolLogP(m), HBD=Lipinski.NumHDonors(m), HBA=Lipinski.NumHAcceptors(m),
                     TPSA=rdMolDescriptors.CalcTPSA(m), RotB=Lipinski.NumRotatableBonds(m), QED=QED.qed(m),
                     heavy=na, Lipinski_viol=int(mw > 500) + int(Descriptors.MolLogP(m) > 5) +
                     int(Lipinski.NumHDonors(m) > 5) + int(Lipinski.NumHAcceptors(m) > 10)))

for n, s in CONTROLS.items(): add(n, s, "control")
for r in df.itertuples():
    nm = (r.name if isinstance(r.name, str) and r.name.strip() else r.chembl_id)
    if nm in lig: nm = f"{nm} ({r.chembl_id})"
    add(nm, r.smiles, r.chembl_id)

props = pd.DataFrame(rows)
props.to_csv(f"{D}/library_properties.csv", index=False)
pickle.dump(lig, open(f"{D}/ligands.pkl", "wb"))
print("prepared ligands:", len(lig))
print(props.describe().loc[["mean", "min", "max"], ["MW", "LogP", "heavy", "RotB"]].round(1))
