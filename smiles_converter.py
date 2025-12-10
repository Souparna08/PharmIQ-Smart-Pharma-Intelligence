import requests
from rdkit import Chem

# Optional quick dictionary for common molecules
NAME_TO_SMILES = {
    "phenol": "c1ccccc1O",
    "benzene": "c1ccccc1",
    "ethanol": "CCO",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O"
}


def fetch_smiles_from_pubchem(name):
    """
    Fetches SMILES string from PubChem using compound name.
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        response = requests.get(url, timeout=5)

        if response.status_code != 200:
            return None

        data = response.json()
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]

    except:
        return None


def convert_name_to_smiles(input_string):
    """
    Detects if the input is name or SMILES.
    Tries predefined dict → PubChem → RDKit validation.
    """
    text = input_string.strip().lower()

    # 1. Check predefined list
    if text in NAME_TO_SMILES:
        return NAME_TO_SMILES[text]

    # 2. Try PubChem API
    pubchem_smiles = fetch_smiles_from_pubchem(text)
    if pubchem_smiles:
        return pubchem_smiles

    # 3. Check if input is already valid SMILES
    mol = Chem.MolFromSmiles(input_string)
    if mol:
        return input_string

    # 4. Completely invalid
    return None
