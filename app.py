import streamlit as st
from patent_agent import load_patent_db, find_similar_patent
from report_generator import generate_pdf_report
from smiles_converter import convert_name_to_smiles


df = load_patent_db("patent_db.csv")

st.title("🔍 PharmIQ - Smart Pharma Intelligence")


smiles_input = st.text_input("Enter SMILES or Molecule Name:")

# Create session state for result
if "result" not in st.session_state:
    st.session_state.result = None


# Check Patent Button

if st.button("Check Patent Status"):
    if not smiles_input.strip():
        st.error("Please enter a valid input.")
    else:
        smiles = convert_name_to_smiles(smiles_input)

        if smiles is None:
            st.error("❌ Invalid molecule name or SMILES string.")
            st.session_state.result = None
        else:
            st.info(f"Detected SMILES: **{smiles}**")

            result = find_similar_patent(smiles, df)
            st.session_state.result = result  # SAVE RESULT

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


# PDF Download Button 

if st.session_state.result is not None:
    if st.button("Generate PDF Report"):
        pdf_file = generate_pdf_report(st.session_state.result, "patent_result.pdf")

        with open(pdf_file, "rb") as f:
            st.download_button(
                label="📄 Download Patent PDF Report",
                data=f.read(),
                file_name="patent_result.pdf",
                mime="application/pdf"
            )
