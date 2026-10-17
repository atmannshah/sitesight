from io import StringIO

import pandas as pd
import requests
from Bio.Data.PDBData import residue_sasa_scales
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from uniprot import phosphosites


def analyse(uid):
    info = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}").json()[0]
    pdb = requests.get(info["pdbUrl"]).text
    structure = PDBParser(QUIET=True).get_structure(uid, StringIO(pdb))
    ShrakeRupley().compute(structure, level="R")
    residues = {r.id[1]: r for r in structure.get_residues()}
    plddt = pd.Series({i: r["CA"].get_bfactor() for i, r in residues.items()}, name="pLDDT")
    sites = pd.DataFrame(
        [
            {
                "pos": pos,
                "res": residues[pos].get_resname(),
                "pLDDT": plddt[pos],
                "RSA": round(residues[pos].sasa / residue_sasa_scales["Sander"][residues[pos].get_resname()], 2),
                "type": kind,
            }
            for pos, kind in phosphosites(uid).items()
        ],
        columns=["pos", "res", "pLDDT", "RSA", "type"],
    )
    return pdb, plddt, sites
