from data_loader.data_loader import DataLoader
from mol_evaluator.evaluator import MolEvaluator


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
    # Load the data
    fake_smiles_path = "../dataset/fake.csv"
    real_smiles_path = "../dataset/real.csv"
    config = {
        "LEVENSHTEIN_THRESHOLD": 0.5,
        "VERY_HIGH_SIMILARITY_THRESHOLD": 0.9,
        "HIGH_SIMILARITY_THRESHOLD": 0.88,
        "LOW_SIMILARITY_THRESHOLD": 0.3,
        "SOLUBILITY_THRESHOLDS": {
            "VERY_HIGH": -1,
            "HIGH": 0,
            "MODERATE": 2,
            "LOW": 4,
            "VERY_LOW": "Infinity"
        },
        "RELEVANT_DESCRIPTORS": [
            "fr_Al_COO", "fr_NH1", "fr_ketone", "fr_halogen",
            "MaxEStateIndex", "MinEStateIndex", "MinPartialCharge", "MaxPartialCharge",
            "fr_COO", "fr_Ar_N", "fr_Ar_OH",
            "MolWt", "ExactMolWt", "HeavyAtomCount", "NumRotatableBonds",
            "FractionCSP3", "LabuteASA", "RingCount",
            "MolLogP", "TPSA",
            "SlogP_VSA1", "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4",
            "SlogP_VSA5", "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "SlogP_VSA9", "SlogP_VSA10",
            "PEOE_VSA1", "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5", "PEOE_VSA6",
            "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12",
            "PEOE_VSA13", "PEOE_VSA14",
            "NumAromaticRings", "NumSaturatedRings", "fr_benzene", "fr_bicyclic",
            "Chi0", "Chi0n", "Chi0v", "Chi1", "Chi1n", "Chi1v",
            "Chi2n", "Chi2v", "Chi3n", "Chi3v", "Chi4n", "Chi4v", "HallKierAlpha"
        ],
        "TANIMOTO_THRESHOLDS": {
            "VERY_HIGH": 0.9,
            "HIGH": 0.88,
            "MODERATE": 0.3
        },
        "VALID_SOLUBILITY_LABELS": ["VERY_HIGH", "HIGH", "MODERATE"],
        "VALID_TANIMOTO_LABELS": ["HIGH", "MODERATE", "LOW"],
        "MAX_SUBSTRUCTURES_MATCHES": 0,
        "REPORT_FOLDER": "./report"
    }

    # Load the data
    data_loader = load_data(real_smiles_path=real_smiles_path, fake_smiles_path=fake_smiles_path)

    # Evaluate the data
    evaluate(data_loader, config)


if __name__ == "__main__":
    main()
