# SiteSight

Fetch an AlphaFold model by UniProt ID and plot per-residue confidence (pLDDT).
[![tests](https://github.com/atmannshah/sitesight/actions/workflows/tests.yml/badge.svg)](https://github.com/atmannshah/sitesight/actions)
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://sitesight.streamlit.app)
## Usage


# SiteSight

Map phosphorylation sites onto AlphaFold structures, with per-site confidence (pLDDT) and surface exposure (RSA).

## Usage

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Why "SiteSight"?

`PhosphositeStructureConfidenceExposureViewer` didn't fit on the logo.

## Caveats

AlphaFold models the unphosphorylated protein. SiteSight shows where the phosphate goes; the phosphate is not included.

```bash
pip install -r requirements.txt
python alphafold.py P04637
```
