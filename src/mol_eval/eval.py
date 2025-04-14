import argparse

from mol_eval.version import __version__
from mol_eval.data_loader import DataLoader
from mol_eval.evaluator import MolEvaluator
from mol_eval.commons import load_config_file
from mol_eval.config import DefaultConfig, ArgsParserConfig
from mol_eval.schemas import ConfigSchema


def evaluate(dl: DataLoader, config: ConfigSchema):
    mol_evaluator = MolEvaluator()
    mol_evaluator.evaluate(dl=dl, config=config)


def parse_args():
    parser = argparse.ArgumentParser(
        description=ArgsParserConfig.tool_description,
    )

    parser.add_argument(
        "--real_data",
        required=False,
        help="Path to the real SMILES data file (CSV).",
        default=DefaultConfig.real_data_path,
    )
    parser.add_argument(
        "--fake_data",
        required=False,
        help="Path to the fake SMILES data file (CSV).",
        default=DefaultConfig.fake_data_path,
    )
    parser.add_argument(
        "--configs",
        required=False,
        help="Path to the configuration JSON file.",
        default=DefaultConfig.config_path,
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    config: ConfigSchema = load_config_file(file_path=args.configs, read_op_type="r")

    data_loader = DataLoader(args.real_data, fake_smiles_path=args.fake_data)
    data_loader.load_smiles()

    mol_evaluator = MolEvaluator()
    mol_evaluator.evaluate(dl=data_loader, config=config)


if __name__ == "__main__":
    main()
