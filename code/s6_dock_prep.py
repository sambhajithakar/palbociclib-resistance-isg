import gemmi, numpy as np, subprocess, json, os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
OUT = "/home/claude/an/dock"; os.makedirs(OUT, exist_ok=True)
PDB = "/mnt/user-data/uploads/Documents/Palbociclib_Resistance_Study/data/raw/1BF5.pdb"
pocket_res = [584, 602, 603, 604, 605, 606, 607, 613, 630, 631, 632, 633]

# crystal pTyr701 of the symmetry mate, expressed in chain A pocket frame (for pose validation)
st = gemmi.read_structure(PDB); m = st[0]
ns = gemmi.NeighborSearch(m, st.cell, 6).populate()
A = m["A"]
pocket_atoms = np.array([[