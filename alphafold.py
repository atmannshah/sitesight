import sys
from io import StringIO

import matplotlib.pyplot as plt
import requests
from Bio.Data.PDBData import residue_sasa_scales
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from uniprot import phosphosites

uid = sys.argv[1]
info = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}").json()[0]
pdb = requests.get(info["pdbUrl"]).text
structure = PDBParser(QUIET=True).get_structure(uid, StringIO(pdb))
ShrakeRupley().compute(structure, level="R")
residues = {r.id[1]: r for r in structure.get_residues()}
plddt = {i: r["CA"].get_bfactor() for i, r in residues.items()}
sites = phosphosites(uid)

print(f"{uid}: {len(plddt)} residues, mean pLDDT {sum(plddt.values()) / len(plddt):.1f}")
print("pos\tres\tpLDDT\tRSA\ttype")
for pos, kind in sites.items():
    r = residues[pos]
    rsa = r.sasa / residue_sasa_scales["Sander"][r.get_resname()]
    print(f"{pos}\t{r.get_resname()}\t{plddt[pos]:.1f}\t{rsa:.2f}\t{kind}")

plt.plot(list(plddt), list(plddt.values()))
plt.scatter(list(sites), [plddt[p] for p in sites], c="red", zorder=3)
plt.axhline(70, ls="--", c="grey")
plt.xlabel("Residue")
plt.ylabel("pLDDT")
plt.title(uid)
plt.savefig(f"{uid}_plddt.png")
