from io import StringIO

import numpy as np
import pandas as pd
import requests
from Bio.Data.PDBData import residue_sasa_scales
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

from uniprot import ptms

CROSSTALK_CUTOFF = 8.0


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
    beta = {r.id[1]: (r["CB"] if "CB" in r else r["CA"]).coord for r in structure.get_residues()}
    sites = pd.DataFrame(ptms(uid), columns=["pos", "ptm", "type"]).merge(residues.reset_index(), on="pos")
    xyz = np.array([beta[p] for p in sites["pos"]]).reshape(-1, 3)
    dist = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    tags = [f"{seq1(r)}{p} {t}" for p, r, t in zip(sites["pos"], sites["res"], sites["ptm"])]
    sites["nearby"] = [
        ", ".join(dict.fromkeys(tag for tag, p, d in zip(tags, sites["pos"], row) if p != pos and d <= CROSSTALK_CUTOFF))
        for pos, row in zip(sites["pos"], dist)
    ]
    return pdb, residues, sites[["pos", "res", "pLDDT", "RSA", "ptm", "type", "nearby"]]
