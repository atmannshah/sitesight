<p align="center"><img src="assets/logo.png" width="180"></p>

<h1 align="center">SiteSight</h1>

<p align="center">
  <a href="https://github.com/atmannshah/sitesight/actions"><img src="https://github.com/atmannshah/sitesight/actions/workflows/tests.yml/badge.svg"></a>
  <a href="https://<your-app>.streamlit.app"><img src="https://static.streamlit.io/badges/streamlit_badge_black_white.svg"></a>
</p>

Map known phosphorylation sites onto AlphaFold structures and check how well each site is modelled and how exposed it is.

## What it does

Enter a UniProt ID and SiteSight:

1. Downloads the predicted structure from [AlphaFold DB](https://alphafold.ebi.ac.uk)
2. Fetches annotated phosphosites from [UniProt](https://www.uniprot.org)
3. Reports, for each site:

| Column | Meaning |
|---|---|
| `pos` | Residue number |
| `res` | Residue (SER, THR, TYR) |
| `pLDDT` | AlphaFold confidence, 0–100. Above 70 is reliable; below 50 is usually disordered |
| `RSA` | Relative solvent accessibility, 0–1. Below 0.25 is buried; above is surface-exposed |
| `type` | Modification, e.g. Phosphoserine |

4. Shows the structure in 3D, coloured by pLDDT, with phosphosites as red sticks

## Run locally

```bash
git clone https://github.com/atmannshah/sitesight.git
cd sitesight
pip install -r requirements.txt
streamlit run app.py
```

## Project structure

```
app.py          Streamlit interface
sitesight.py    Structure download, pLDDT and RSA calculation
uniprot.py      Phosphosite lookup
tests/          Unit tests (APIs mocked)
```

## Tests

```bash
python -m pytest
```

## Limitations

- AlphaFold models the unmodified protein, so phosphate groups are not present in the structure.
- RSA is computed on the isolated monomer; sites buried in complexes may appear exposed.
- Only sites annotated in UniProt are shown.

## License

MIT
