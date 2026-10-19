import streamlit as st

st.set_page_config(page_title="SiteSight", layout="wide")
st.logo("assets/wordmark.svg", size="large")
st.navigation([st.Page("explore.py", title="Explore", default=True), st.Page("about.py", title="About")], position="top").run()
