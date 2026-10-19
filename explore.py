import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from sitesight import analyse

st.caption("Phosphates can't run, they can't hide either...")
uid = st.text_input("UniProt ID", "P04637", label_visibility="collapsed", placeholder="UniProt ID, e.g. P04637").strip().upper()

try:
    with st.spinner("Pinpointing those phosphosites for you..."):
        pdb, plddt, sites = st.cache_data(analyse)(uid)
except Exception:
    st.error(f"No AlphaFold model for {uid}. Either it's a typo, or it's still folding.")
    st.stop()

a, b, c = st.columns(3)
a.metric("Residues", len(plddt))
b.metric("Mean pLDDT", f"{plddt.mean():.1f}")
c.metric("Phosphosites", len(sites))

left, right = st.columns([3, 2], gap="large")

with left:
    view = py3Dmol.view(width=700, height=520)
    view.addModel(pdb, "pdb")
    view.setBackgroundColor("white")
    view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}})
    if len(sites):
        view.addStyle({"resi": sites["pos"].tolist()}, {"stick": {"colorscheme": "redCarbon"}})
    view.zoomTo()
    components.html(view._make_html(), height=520)
    st.caption("Spaghetti regions are low pLDDT. Not a bad prediction, just a protein with commitment issues.")

with right:
    st.line_chart(plddt, height=220)
    if sites.empty:
        st.info("No annotated phosphosites. The kinases have left this one on read.")
    else:
        st.dataframe(sites, hide_index=True)

st.caption("Structures from AlphaFold DB, sites from UniProt. SiteSight just points at things.")
