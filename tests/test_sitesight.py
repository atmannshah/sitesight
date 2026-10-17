from unittest.mock import Mock

import sitesight
import uniprot

PDB = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C
ATOM      2  CA  SER A   2       3.800   0.000   0.000  1.00 60.00           C
ATOM      3  CA  GLY A   3       7.600   0.000   0.000  1.00 30.00           C
END
"""

FEATURES = {
    "features": [
        {"type": "Modified residue", "location": {"start": {"value": 2}}, "description": "Phosphoserine; by ATM"},
        {"type": "Modified residue", "location": {"start": {"value": 3}}, "description": "N6-acetyllysine"},
        {"type": "Helix", "location": {"start": {"value": 1}}, "description": ""},
    ]
}


def fake_get(url, *args, **kwargs):
    if "api/prediction" in url:
        return Mock(json=lambda: [{"pdbUrl": "model.pdb"}])
    if "uniprot" in url:
        return Mock(json=lambda: FEATURES)
    return Mock(text=PDB)


def test_only_phosphates_get_an_invite(monkeypatch):
    monkeypatch.setattr(uniprot.requests, "get", fake_get)
    assert uniprot.phosphosites("X") == {2: "Phosphoserine"}


def test_analyse(monkeypatch):
    monkeypatch.setattr(sitesight.requests, "get", fake_get)
    monkeypatch.setattr(uniprot.requests, "get", fake_get)
    pdb, plddt, sites = sitesight.analyse("X")
    assert plddt.tolist() == [90, 60, 30]
    assert sites[["pos", "res", "pLDDT"]].values.tolist() == [[2, "SER", 60]]
    assert 0 < sites["RSA"][0]
