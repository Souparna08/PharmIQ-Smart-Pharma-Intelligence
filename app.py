import streamlit as st
import pandas as pd
from patent_agent import load_patent_db, find_similar_patent
from report_generator import generate_pdf_report
from smiles_converter import convert_name_to_smiles   # ← NEW

# Load DB
df = load_patent_db("patent_db.csv")

st.title("🔍 PharmIQ - Smart Pharma Intelligence")
st.write("Enter a molecule SMILES or chemical name (e.g., phenol, aspirin, ethanol).")

smiles_input = st.text_input("Enter SMILES or Molecule Name:")

if st.button("Check Patent Status"):
    if smiles_input.strip() == "":
        st.error("Please enter a valid input.")
    else:
        # Convert name → SMILES
        smiles = convert_name_to_smiles(smiles_input)

    if smiles_input.lower() != smiles:
        st.info(f"Detected SMILES for **{smiles_input}** → `{smiles}`")

        if smiles is None:
            st.error("❌ Invalid molecule name or SMILES string.")
        else:
            st.info(f"Detected SMILES: **{smiles}**")

            # Run similarity
            result = find_similar_patent(smiles, df)

            if result is None:
                st.success("✅ No similar patent found! Molecule seems novel.")
            else:
                st.warning("⚠️ Similar patented molecule found!")
                st.write("### 🔬 Closest Patent Match")
                st.write(f"**Patent ID:** {result['patent_id']}")
                st.write(f"**Title:** {result['title']}")
                st.write(f"**Year:** {result['year']}")
                st.write(f"**Similarity Score:** {result['similarity']:.2f}")
                st.write(f"**URL:** {result['url']}")
                st.write("---")

                if st.button("Generate PDF Report"):
                    pdf_file = generate_pdf_report(result, "patent_result.pdf")
                    with open(pdf_file, "rb") as f:
                        st.download_button(
                            label="📄 Download Patent PDF Report",
                            data=f,
                            file_name="patent_result.pdf",
                            mime="application/pdf"
                        )
