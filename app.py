import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from sitesight import analyse

st.set_page_config(page_title="SiteSight", layout="wide")
st.title("SiteSight")
st.caption("Kinases know exactly where they're going. Now you can too.")
uid = st.text_input("UniProt ID", "P04637").strip().upper()

try:
    with st.spinner("Folding... kidding, AlphaFold already did that."):
        pdb, plddt, sites = st.cache_data(analyse)(uid)
except Exception:
    st.error(f"No AlphaFold model for {uid}. Either it's a typo, or it's still folding.")
    st.stop()

left, right = st.columns(2)

with left:
    view = py3Dmol.view(width=600, height=500)
    view.addModel(pdb, "pdb")
    view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}})
    if len(sites):
        view.addStyle({"resi": sites["pos"].tolist()}, {"stick": {"colorscheme": "redCarbon"}})
    view.zoomTo()
    components.html(view._make_html(), height=500)
    st.caption("Spaghetti regions are low pLDDT. Not a bad prediction, just a protein with commitment issues.")

with right:
    st.metric("Mean pLDDT", f"{plddt.mean():.1f}")
    st.line_chart(plddt)
    if sites.empty:
        st.info("No annotated phosphosites. The kinases have left this one on read.")
    else:
        st.dataframe(sites, hide_index=True)

st.caption("Structures from AlphaFold DB, sites from UniProt. SiteSight just points at things.")
