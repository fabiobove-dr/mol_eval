from dataclasses import dataclass


@dataclass
class DefaultConfig:
    real_data_path: str = "../dataset/real.csv"
    fake_data_path: str = "../dataset/fake.csv"
    config_path: str = "../dataset/config.json"


@dataclass
class ArgsParserConfig:
    tool_description: str = "Molecule Evaluator: Evaluate real and fake SMILES data using a configuration file."
