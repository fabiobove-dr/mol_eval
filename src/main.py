import json

from rdkit import Chem

from src.data_loader.data_loader import DataLoader
from src.mol_evaluator.mol_evaluator import MolEvaluator


def load_config(file_path):
    with open(file_path, 'r') as file:
        _config = json.load(file)
    return _config


# Load the configuration
config = load_config('../examples/configs.json')

# Access configuration variables
LEVENSHTEIN_THRESHOLD = config["LEVENSHTEIN_THRESHOLD"]
VERY_HIGH_SIMILARITY_THRESHOLD = config["VERY_HIGH_SIMILARITY_THRESHOLD"]
HIGH_SIMILARITY_THRESHOLD = config["HIGH_SIMILARITY_THRESHOLD"]
LOW_SIMILARITY_THRESHOLD = config["LOW_SIMILARITY_THRESHOLD"]
SOLUBILITY_THRESHOLDS = config["SOLUBILITY_THRESHOLDS"]
RELEVANT_DESCRIPTORS = config["RELEVANT_DESCRIPTORS"]
TANIMOTO_THRESHOLDS = config["TANIMOTO_THRESHOLDS"]
VALID_SOLUBILITY_LABELS = config["VALID_SOLUBILITY_LABELS"]
VALID_TANIMOTO_LABELS = config["VALID_TANIMOTO_LABELS"]
MAX_SUBSTRUCTURES_MATCHES = config["MAX_SUBSTRUCTURES_MATCHES"]
REPORT_FOLDER = config["REPORT_FOLDER"]


def load_data(real_smiles_path: str, fake_smiles_path: str):
    dl = DataLoader(real_smiles_path, fake_smiles_path)
    dl.load_smiles()
    print("Data loaded successfully!")
    print(f"Real SMILES Size: {len(dl.get_real_smiles())}")
    print(f"Fake SMILES Size: {len(dl.get_fake_smiles())}")
    return dl


def evaluate(dl: DataLoader):
    mol_evaluator = MolEvaluator()

    print("Removing non-molecules...")
    df = mol_evaluator.remove_non_molecules(dl.fake_smiles_df)
    print(f"Valid molecules: {len(df)}")

    print("Removing existing smiles...")
    df = mol_evaluator.remove_existing(df, dl.real_smiles_df)
    print(f"Valid molecules: {len(df)}")

    print("Removing duplicate smiles...")
    df = mol_evaluator.remove_duplicates(df)
    print(f"Valid molecules: {len(df)}")

    print("Adding similarity with Levenshtein...")
    df = mol_evaluator.add_levenshtein_similarity(df, dl.real_smiles_df, threshold=LEVENSHTEIN_THRESHOLD)
    print(f"Valid molecules: {len(df)}")

    print("Adding Descriptors...")
    df = mol_evaluator.describe_fake_smiles(df, RELEVANT_DESCRIPTORS)
    print(f"Valid molecules: {len(df)}")

    print("Adding Water Solubility Label...")
    df = mol_evaluator.add_solubility_labels(df, SOLUBILITY_THRESHOLDS)
    print(f"Valid molecules: {len(df)}")

    print("Filtering molecules with low solubility...")
    df = mol_evaluator.filter_by_solubility(df, VALID_SOLUBILITY_LABELS)
    print(f"Valid molecules: {len(df)}")

    print("Computing substructure matches...")
    df = mol_evaluator.compute_substructure_matches(df, dl.real_smiles_df)
    print(f"Valid molecules: {len(df)}")

    print("Filtering molecules with high substructure matches...")
    df = mol_evaluator.filter_by_substructure_matches_number(df, MAX_SUBSTRUCTURES_MATCHES)
    print(f"Valid molecules: {len(df)}")

    print("Computing Tanimoto Similarity...")
    df = mol_evaluator.add_tanimoto_similarity_score_and_label(df, dl.real_smiles_df, TANIMOTO_THRESHOLDS)
    print(f"Valid molecules: {len(df)}")

    print("Filtering molecules with high Tanimoto Similarity...")
    df = mol_evaluator.filter_by_tanimoto_label(df, VALID_TANIMOTO_LABELS)
    print(f"Valid molecules: {len(df)}")

    print("Adding Image...")
    df = mol_evaluator.add_2d_visualizations(df)
    print(f"Valid molecules: {len(df)}")

    # TODO: move this to eval and refactor.
    print("Adding PatCID Labels...")
    PATCHID_PATH = "../dataset/patcid_molecule_to_patents_fixed.jsonl"
    all_smiles = df["smiles"].values.tolist()
    # Open the file and search for the SMILES query
    founds = []
    patents = []

    with open(PATCHID_PATH, "r", encoding="utf-8") as file:
        for smiles in all_smiles:
            query_smiles = smiles
            # Canonicalize SMILES
            molecule = Chem.MolFromSmiles(query_smiles)
            query_smiles = Chem.MolToSmiles(molecule)
            query_mol = Chem.MolFromSmiles(query_smiles)
            query_inchikey = Chem.MolToInchiKey(query_mol)

            # Search SMILES in PatCID (Based on this naive approach a search takes order-of-magnitude 10 seconds.)
            # Initialize variables for result and error handling
            smiles_entry = None
            file.seek(0)
            for line in file:
                if query_inchikey in line or query_smiles in line:
                    smiles_entry = json.loads(line)  # Assuming the file contains valid JSON on each line
                    break

            # Check if SMILES was found
            if smiles_entry is None:
                founds.append("PATENT_NOT_FOUND")
                patents.append([])
            else:
                founds.append("PATENT_FOUND")
                patents.append([p["id"] for p in smiles_entry["patents"]])

        df["patents_labels"] = founds
        df["patents"] = patents

    print(f"Saving {df.shape[0]} results...")
    mol_evaluator.create_report(df, REPORT_FOLDER)
    print("Done!")


if __name__ == "__main__":
    data_loader = load_data(real_smiles_path="../dataset/original.csv", fake_smiles_path="../dataset/fake.csv")
    evaluate(data_loader)
