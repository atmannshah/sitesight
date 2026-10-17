import requests


def phosphosites(uid):
    entry = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.json").json()
    return {
        f["location"]["start"]["value"]: f["description"].split(";")[0]
        for f in entry.get("features", [])
        if f["type"] == "Modified residue" and f["description"].startswith("Phospho")
    }
