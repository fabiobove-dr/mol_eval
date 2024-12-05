from pymongo import MongoClient
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

# MongoDB connection details
username = "root"
password = "example"
host = "localhost"
port = 27017
db_name = "moleculeDB"
auth_db = "admin"  # Admin database for authentication

# MongoDB connection with authentication
client = MongoClient(f"mongodb://{username}:{password}@{host}:{port}/?authSource={auth_db}")
db = client[db_name]
collection = db['patcid']  # You can change this to the desired collection

# Load the list of SMILES from a CSV file
df = pd.read_csv('../report/report.csv')
list_of_smiles = df['smiles'].tolist()

# Initialize a dictionary to store SMILES and their associated patents
smiles_to_patents = {}

# Iterate over the list of SMILES
for smiles in list_of_smiles:
    try:
        # Canonicalize SMILES
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        canonical_smiles = Chem.MolToSmiles(molecule)
        query_mol = Chem.MolFromSmiles(canonical_smiles)
        query_inchikey = Chem.MolToInchiKey(query_mol)

        # Query the database by InChIKey or SMILES
        results = collection.find({
            "$or": [
                {"molecule.smiles": {"$regex": canonical_smiles, "$options": "i"}},
                {"molecule.inchikey": query_inchikey}
            ]
        })

        # Collect the list of patent IDs
        patent_ids = []
        for result in results:
            patents = result.get('patents', [])  # Default to empty if 'patents' field is missing
            if isinstance(patents, list):  # If 'patents' is a list of dictionaries
                patent_ids.extend(patent.get('id') for patent in patents if 'id' in patent)
            elif isinstance(patents, dict):  # If 'patents' is a single dictionary
                patent_ids.append(patents.get('id'))

        # Add to the dictionary (use an empty list if no patents found)
        smiles_to_patents[smiles] = patent_ids if patent_ids else []

    except Exception as e:
        print(f"Error processing SMILES '{smiles}': {e}")
        smiles_to_patents[smiles] = []  # Add empty list in case of an error

output_df = pd.DataFrame({
    'SMILES': list(smiles_to_patents.keys()),
    'Patent IDs': list(smiles_to_patents.values())
})

# Save to CSV
output_df.to_csv('../report/smiles_to_patents.csv', index=False)
