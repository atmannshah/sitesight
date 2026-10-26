import json

import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from Bio.SeqUtils import seq1

from sitesight import analyse

CLICK = """function(atom, viewer) {
    var rsa = RSA;
    viewer.picked = viewer.picked || {};
    if (viewer.picked[atom.resi]) {
        viewer.removeLabel(viewer.picked[atom.resi]);
        delete viewer.picked[atom.resi];
    } else {
        viewer.picked[atom.resi] = viewer.addLabel(
            atom.resn + " " + atom.resi + "  pLDDT " + atom.b.toFixed(1) + "  RSA " + rsa[atom.resi],
            {position: atom, fontSize: 12, fontColor: "white", backgroundColor: "#0053D6", backgroundOpacity: 0.9}
        );
    }
    viewer.render();
}"""

st.caption("Kinases know exactly where they're going. Now you can too.")
uid = st.text_input("UniProt ID", "P04637", label_visibility="collapsed", placeholder="UniProt ID, e.g. P04637").strip().upper()

try:
    with st.spinner("Folding... kidding, AlphaFold already did that."):
        pdb, residues, sites = st.cache_data(analyse)(uid)
except Exception:
    st.error(f"No AlphaFold model for {uid}. Either it's a typo, or it's still folding.")
    st.stop()

a, b, c = st.columns(3)
a.metric("Residues", len(residues))
b.metric("Mean pLDDT", f"{residues['pLDDT'].mean():.1f}")
c.metric("PTM sites", len(sites))

left, right = st.columns([3, 2], gap="large")

with right:
    st.line_chart(residues["pLDDT"], height=220)
    if sites.empty:
        st.info("No annotated PTMs. This protein is exactly as translated, apparently.")
        picked = None
    else:
        rows = st.dataframe(sites, hide_index=True, on_select="rerun", selection_mode="single-row").selection.rows
        picked = int(sites["pos"][rows[0]]) if rows else None
        st.caption("Select a row to zoom to that site.")

with left:
    mode = st.segmented_control("View", ["Confidence", "PTMs"], default="Confidence", label_visibility="collapsed") or "Confidence"
    view = py3Dmol.view(width=700, height=520)
    view.addModel(pdb, "pdb")
    view.setBackgroundColor("white")
    if mode == "Confidence":
        view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}})
        caption = "Blue is confident (pLDDT > 90), red is not (< 50). Spaghetti regions are a protein with commitment issues."
    else:
        view.setStyle({"cartoon": {"color": "#D1D5DB"}})
        for pos, res in zip(sites["pos"].tolist(), sites["res"]):
            view.addStyle({"resi": pos}, {"stick": {"colorscheme": "redCarbon", "radius": 0.3}})
            view.addLabel(f"{seq1(res)}{pos}", {"fontSize": 11, "fontColor": "#111827", "backgroundColor": "white", "backgroundOpacity": 0.85, "borderThickness": 0}, {"resi": pos, "atom": "CA"})
        caption = "PTM sites shown as red sticks, labelled by residue and position."
    view.setClickable({}, True, CLICK.replace("RSA;", json.dumps(residues["RSA"].to_dict()) + ";"))
    if picked:
        view.addStyle({"resi": picked}, {"stick": {"colorscheme": "yellowCarbon", "radius": 0.35}})
        view.zoomTo({"within": {"distance": 12, "sel": {"resi": picked}}})
    else:
        view.zoomTo()
    components.html(view._make_html(), height=520)
    st.caption(caption + " Click any residue for details, click again to hide.")

st.caption("Structures from AlphaFold DB, sites from UniProt. SiteSight just points at things.")
