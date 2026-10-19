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
    residues = pd.DataFrame(
        [
            {
                "pos": r.id[1],
                "res": r.get_resname(),
                "pLDDT": r["CA"].get_bfactor(),
                "RSA": round(r.sasa / residue_sasa_scales["Sander"][r.get_resname()], 2),
            }
            for r in structure.get_residues()
        ]
    ).set_index("pos")
    sites = residues.join(pd.Series(phosphosites(uid), name="type", dtype=str), how="inner").rename_axis("pos").reset_index()
    return pdb, residues, sites
