from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd

# ---------------------------
# Load Patent Database
# ---------------------------
def load_patent_db(path):
    return pd.read_csv(path)

# ---------------------------
# Compute Molecular Fingerprint
# ---------------------------
def get_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)

# ---------------------------
# Compare Query Molecule With Patent DB
# ---------------------------
def find_similar_patent(query_smiles, df):
    query_fp = get_fingerprint(query_smiles)

    if query_fp is None:
        return None

    best_score = 0
    best_match = None

    for _, row in df.iterrows():
        patent_fp = get_fingerprint(row["smiles"])
        if patent_fp is None:
            continue
        
        score = DataStructs.TanimotoSimilarity(query_fp, patent_fp)

        if score > best_score:
            best_score = score
            best_match = row

    if best_score < 0.4:     # threshold
        return None

    result = best_match.to_dict()
    result["similarity"] = best_score

    return result
