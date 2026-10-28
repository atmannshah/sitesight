from unittest.mock import Mock

import sitesight
import uniprot

PDB = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C
ATOM      2  CA  SER A   2       3.800   0.000   0.000  1.00 60.00           C
ATOM      3  CA  LYS A   3       7.600   0.000   0.000  1.00 30.00           C
END
"""


def feature(kind, pos, description):
    return {"type": kind, "location": {"start": {"value": pos}}, "description": description}


FEATURES = {
    "features": [
        feature("Modified residue", 2, "Phosphoserine; by ATM"),
        feature("Modified residue", 3, "N6-acetyllysine; alternate"),
        feature("Modified residue", 3, "N6,N6,N6-trimethyllysine; alternate"),
        feature("Cross-link", 3, "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)"),
        feature("Cross-link", 3, "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)"),
        feature("Glycosylation", 2, "O-linked (GlcNAc) serine"),
        feature("Modified residue", 3, "N6-(pyridoxal phosphate)lysine"),
        feature("Helix", 1, ""),
    ]
}


def fake_get(url, *args, **kwargs):
    if "api/prediction" in url:
        return Mock(json=lambda: [{"pdbUrl": "model.pdb"}])
    if "uniprot" in url:
        return Mock(json=lambda: FEATURES)
    return Mock(text=PDB)


def test_every_ptm_gets_the_right_label(monkeypatch):
    monkeypatch.setattr(uniprot.requests, "get", fake_get)
    assert [(p["pos"], p["ptm"]) for p in uniprot.ptms("X")] == [
        (2, "Phospho"),
        (3, "Acetyl"),
        (3, "Methyl"),
        (3, "Ubiquitin"),
        (3, "SUMO"),
        (2, "Glyco"),
        (3, "Other"),
    ]


def test_analyse(monkeypatch):
    monkeypatch.setattr(sitesight.requests, "get", fake_get)
    pdb, residues, sites = sitesight.analyse("X")
    assert residues["pLDDT"].tolist() == [90, 60, 30]
    assert (residues["RSA"] > 0).all()
    assert sites[["pos", "res", "pLDDT", "ptm"]].values.tolist()[:2] == [[2, "SER", 60, "Phospho"], [3, "LYS", 30, "Acetyl"]]
    assert len(sites) == 7


def test_no_sites(monkeypatch):
    monkeypatch.setitem(FEATURES, "features", [])
    monkeypatch.setattr(sitesight.requests, "get", fake_get)
    assert sitesight.analyse("X")[2].empty





def test_crosstalk(monkeypatch):
    monkeypatch.setattr(sitesight.requests, "get", fake_get)
    monkeypatch.setitem(
        FEATURES,
        "features",
        [
            feature("Modified residue", 1, "N-acetylalanine"),
            feature("Modified residue", 2, "Phosphoserine"),
            feature("Modified residue", 3, "N6-acetyllysine"),
            feature("Modified residue", 3, "N6-methyllysine"),
        ],
    )
    monkeypatch.setattr(sitesight, "CROSSTALK_CUTOFF", 5.0)
    nearby = sitesight.analyse("X")[2].set_index("type")["nearby"]
    assert nearby["N-acetylalanine"] == "S2 Phospho"
    assert nearby["Phosphoserine"] == "A1 Acetyl, K3 Acetyl, K3 Methyl"
    assert nearby["N6-acetyllysine"] == "S2 Phospho"    
