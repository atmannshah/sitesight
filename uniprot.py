import requests

KINDS = {"Modified residue", "Cross-link", "Glycosylation"}


def category(feature):
    d = feature["description"].lower()
    if feature["type"] == "Glycosylation":
        return "Glyco"
    if d.startswith("phospho"):
        return "Phospho"
    return next((name for name in ("Acetyl", "Methyl", "Ubiquitin", "SUMO") if name.lower() in d), "Other")


def ptms(uid):
    entry = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json").json()
    return [
        {"pos": f["location"]["start"]["value"], "ptm": category(f), "type": f["description"].split(";")[0]}
        for f in entry.get("features", [])
        if f["type"] in KINDS
    ]
