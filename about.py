import streamlit as st

st.markdown("""
## About

Phosphorylation is one of the most common ways a cell switches a protein on or off. A kinase adds a phosphate group to a serine, threonine or tyrosine, and that small charged group changes which partners the protein binds, where it goes in the cell and how long it survives.

Thousands of phosphosites are catalogued in databases like UniProt, but a residue number on its own says very little. To reason about a site you want to know where it sits in the folded protein. Is it on the surface where a kinase or a binding domain can reach it, or buried in the core? Is that part of the structure well defined, or is it a flexible tail?

SiteSight answers those questions in a few seconds. It places every annotated phosphosite on the AlphaFold model of the protein and reports two numbers for each one.

**pLDDT** is AlphaFold's own confidence in the local structure, from 0 to 100. Above 70 the backbone is usually reliable. Below 50 the region is often intrinsically disordered, which is common for phosphosites since many sit in flexible linkers and tails.

**RSA** (relative solvent accessibility) is the fraction of a residue's surface exposed to solvent, from 0 to 1. Values below about 0.25 mean the residue is buried.

### Background
I like proteins, I like folding, and I like looking for super specific things. So I started working on SiteSight to make my life less tedious, hopefully yours too :-)
SiteSight grew out of my interest in the structural characterisation and de novo design of protein phosphorylation sites. A recurring first step in that work was simply looking at where a phosphosite sits in a structure. This tool makes that step quick and shareable.

### Data

Structures come from the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk) and phosphosite annotations from [UniProt](https://www.uniprot.org). AlphaFold models the unmodified protein, so phosphate groups are not present in the structures shown.

### Source

[github.com/atmannshah/sitesight](https://github.com/atmannshah/sitesight)

###
"Following the foot-steps of a rag doll dance,
 We are entranced, spellbound, spellbound, spellbound"

-Siouxsie and the Banshees
###



""")
