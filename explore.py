import json

import altair as alt
import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from Bio.SeqUtils import seq1

from sitesight import analyse

PTM_COLOURS = {
    "Phospho": ("#DC2626", "🔴"),
    "Acetyl": ("#EA580C", "🟠"),
    "Methyl": ("#16A34A", "🟢"),
    "Ubiquitin": ("#2563EB", "🔵"),
    "SUMO": ("#9333EA", "🟣"),
    "Glyco": ("#92400E", "🟤"),
    "Other": ("#6B7280", "⚫"),
}

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

st.caption("Enzymes know exactly where they're going. Now you can too.")
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

present = [p for p in PTM_COLOURS if p in set(sites["ptm"])]
chosen = present and st.pills("PTM types", present, selection_mode="multi", default=present, format_func=lambda p: f"{PTM_COLOURS[p][1]} {p}", label_visibility="collapsed")
shown = sites[sites["ptm"].isin(chosen)].reset_index(drop=True)

left, right = st.columns([3, 2], gap="large")

with right:
    confidence, landscape = st.tabs(["pLDDT", "PTM landscape"])
    confidence.line_chart(residues["pLDDT"], height=480, x_label="Residue", y_label="pLDDT")
    points = alt.Chart(shown).mark_circle(size=90, opacity=0.8).encode(
        x=alt.X("pLDDT", scale=alt.Scale(domain=[0, 100])),
        y=alt.Y("RSA", scale=alt.Scale(domain=[0, 1], clamp=True)),
        color=alt.Color("ptm", scale=alt.Scale(domain=list(PTM_COLOURS), range=[c for c, _ in PTM_COLOURS.values()]), legend=None),
        tooltip=["pos", "res", "ptm", "type", "pLDDT", "RSA"],
    )
    rules = alt.Chart(pd.DataFrame({"pLDDT": [70], "RSA": [0.25]}))
    chart = points + rules.mark_rule(strokeDash=[4, 4], color="#9CA3AF").encode(x="pLDDT") + rules.mark_rule(strokeDash=[4, 4], color="#9CA3AF").encode(y="RSA")
    landscape.altair_chart(chart.properties(height=480))
    landscape.caption("Dashed lines mark pLDDT 70 and RSA 0.25. Top left is exposed and disordered, which is where most PTMs like to hang out.")

picked = None
if sites.empty:
    st.info("No annotated PTMs. This protein is exactly as translated, apparently.")
elif shown.empty:
    st.info("Pick a PTM type above.")
else:
    rows = st.dataframe(shown, hide_index=True, on_select="rerun", selection_mode="single-row").selection.rows
    picked = int(shown["pos"][rows[0]]) if rows else None
    st.download_button("Download CSV", shown.to_csv(index=False), f"{uid}_ptms.csv", "text/csv")
    if rows and shown["nearby"][rows[0]]:
        st.caption(f"Within 8 Å of {seq1(shown['res'][rows[0]])}{picked}: {shown['nearby'][rows[0]]}. Possible crosstalk, or just good neighbours.")
    else:
        st.caption("Select a row to zoom to that site. The nearby column lists other PTMs within 8 Å.")

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
        for pos, ptm in zip(shown["pos"].tolist(), shown["ptm"]):
            view.addStyle({"resi": pos}, {"stick": {"color": PTM_COLOURS[ptm][0], "radius": 0.3}})
        for pos, res in shown.drop_duplicates("pos")[["pos", "res"]].values.tolist():
            view.addLabel(f"{seq1(res)}{pos}", {"fontSize": 11, "fontColor": "#111827", "backgroundColor": "white", "backgroundOpacity": 0.85, "borderThickness": 0}, {"resi": pos, "atom": "CA"})
        caption = "PTM sites shown as sticks, coloured by type and labelled by residue and position."
    view.setClickable({}, True, CLICK.replace("RSA;", json.dumps(residues["RSA"].to_dict()) + ";"))
    if picked:
        view.addStyle({"resi": picked}, {"stick": {"colorscheme": "yellowCarbon", "radius": 0.35}})
        view.zoomTo({"within": {"distance": 12, "sel": {"resi": picked}}})
    else:
        view.zoomTo()
    components.html(view._make_html(), height=520)
    st.caption(caption + " Click any residue for details, click again to hide.")

st.caption("Structures from AlphaFold DB, sites from UniProt. SiteSight just points at things.")
