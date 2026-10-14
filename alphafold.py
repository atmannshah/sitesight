import sys
from io import StringIO

import matplotlib.pyplot as plt
import requests
from Bio.PDB import PDBParser

uid = sys.argv[1]
info = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}").json()[0]
pdb = requests.get(info["pdbUrl"]).text
structure = PDBParser(QUIET=True).get_structure(uid, StringIO(pdb))
plddt = [r["CA"].get_bfactor() for r in structure.get_residues()]

print(f"{uid}: {len(plddt)} residues, mean pLDDT {sum(plddt) / len(plddt):.1f}")

plt.plot(range(1, len(plddt) + 1), plddt)
plt.axhline(70, ls="--", c="grey")
plt.xlabel("Residue")
plt.ylabel("pLDDT")
plt.title(uid)
plt.savefig(f"{uid}_plddt.png")
