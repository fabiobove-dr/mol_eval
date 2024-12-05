import argparse
import json

from data_loader.data_loader import DataLoader
from mol_evaluator.mol_evaluator import MolEvaluator


def load_config(file_path):
    """Load configuration from a JSON file."""
    with open(file_path, 'r') as file:
        return json.load(file)


def load_data(real_smiles_path: str, fake_smiles_path: str):
    """Load real and fake SMILES data."""
    dl = DataLoader(real_smiles_path, fake_smiles_path)
    dl.load_smiles()
    print("Data loaded successfully!")
    print(f"Real SMILES Size: {len(dl.get_real_smiles())}")
    print(f"Fake SMILES Size: {len(dl.get_fake_smiles())}")
    return dl


def evaluate(dl: DataLoader, config):
    """Evaluate molecules based on the configuration."""
    mol_evaluator = MolEvaluator()

    thresholds = {
        "levenshtein": config["LEVENSHTEIN_THRESHOLD"],
        "tanimoto": config["TANIMOTO_THRESHOLDS"],
        "solubility": config["SOLUBILITY_THRESHOLDS"],
        "valid_solubility": config["VALID_SOLUBILITY_LABELS"],
        "valid_tanimoto": config["VALID_TANIMOTO_LABELS"],
        "max_substructures": config["MAX_SUBSTRUCTURES_MATCHES"],
    }

    descriptors = config["RELEVANT_DESCRIPTORS"]
    report_folder = config["REPORT_FOLDER"]

    print("🧑‍🔬 Removing non-molecules...")
    df = mol_evaluator.remove_non_molecules(dl.fake_smiles_df)
    print(f"🔬 Valid molecules: {len(df)}")

    print("💀 Removing existing smiles...")
    df = mol_evaluator.remove_existing(df, dl.real_smiles_df)
    print(f"✔️ Valid molecules: {len(df)}")

    print("🔄 Removing duplicate smiles...")
    df = mol_evaluator.remove_duplicates(df)
    print(f"✅ Valid molecules: {len(df)}")

    print("🔍 Adding similarity with Levenshtein...")
    df = mol_evaluator.add_levenshtein_similarity(df, dl.real_smiles_df, threshold=thresholds["levenshtein"])
    print(f"🧑‍🔬 Valid molecules: {len(df)}")

    print("📊 Adding Descriptors...")
    df = mol_evaluator.describe_fake_smiles(df, descriptors)
    print(f"🧬 Valid molecules: {len(df)}")

    print("💧 Adding Water Solubility Label...")
    df = mol_evaluator.add_solubility_labels(df, thresholds["solubility"])
    print(f"🌊 Valid molecules: {len(df)}")

    print("🚫 Filtering molecules with low solubility...")
    df = mol_evaluator.filter_by_solubility(df, thresholds["valid_solubility"])
    print(f"🌱 Valid molecules: {len(df)}")

    print("🧩 Computing substructure matches...")
    df = mol_evaluator.compute_substructure_matches(df, dl.real_smiles_df)
    print(f"🔍 Valid molecules: {len(df)}")

    print("⚖️ Filtering molecules with high substructure matches...")
    df = mol_evaluator.filter_by_substructure_matches_number(df, thresholds["max_substructures"])
    print(f"🔑 Valid molecules: {len(df)}")

    print("📏 Computing Tanimoto Similarity...")
    df = mol_evaluator.add_tanimoto_similarity_score_and_label(df, dl.real_smiles_df, thresholds["tanimoto"])
    print(f"💡 Valid molecules: {len(df)}")

    print("🚫 Filtering molecules with high Tanimoto Similarity...")
    df = mol_evaluator.filter_by_tanimoto_label(df, thresholds["valid_tanimoto"])
    print(f"🧪 Valid molecules: {len(df)}")

    print("🎨 Adding Image...")
    df = mol_evaluator.add_2d_visualizations(df)
    print(f"🖼️ Valid molecules: {len(df)}")

    print(f"🛟 Saving {df.shape[0]} results...")
    mol_evaluator.create_report(df, report_folder)
    print("🥳 Done!")


def main():
    parser = argparse.ArgumentParser(
        description="Molecule Evaluator: Evaluate real and fake SMILES data using a configuration file."
    )
    parser.add_argument("--real_data", required=True, help="Path to the real SMILES data file (CSV).")
    parser.add_argument("--fake_data", required=True, help="Path to the fake SMILES data file (CSV).")
    parser.add_argument("--configs", required=True, help="Path to the configuration JSON file.")

    args = parser.parse_args()

    # Load the configuration file
    config = load_config(args.configs)

    # Load the data
    data_loader = load_data(real_smiles_path=args.real_data, fake_smiles_path=args.fake_data)

    # Evaluate the data
    evaluate(data_loader, config)


if __name__ == "__main__":
    main()
