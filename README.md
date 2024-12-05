<figure>
  <img src="icon.png" alt="description" >
</figure>

## Molecules Evaluator: a tool for the evaluation of molecules smiles. a.k.a. `mol_eval`

<hr>

[![Coverage](https://codecov.io/github/tacclab/bio_dataset_manager/coverage.svg?branch=main)](https://codecov.io/gh/tacclab/bio_dataset_manager) 
[![PyPI Latest Release](https://img.shields.io/pypi/v/bio_dataset_manager.svg)](https://pypi.org/project/bio_dataset_manager/)
![Unit Tests](https://github.com/tacclab/bio_dataset_manager/actions/workflows/main.yml/badge.svg)<br>
[![Powered by TaccLab](https://img.shields.io/badge/powered%20by-TaccLab-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://tacclab.org)<br> 
[![License](https://img.shields.io/github/license/tacclab/bio_dataset_manager.svg)](https://github.com/tacclab/bio_dataset_manager/blob/main/LICENSE)<br>


## Configuration

1. Prepare your dataset files:
- `original.csv`: Must contain two columns `cmpd_name`, `smile`. The columns contain the information on smiles names and sequences.
- `fake.csv`: Must contain one column `smile`, that contains the smiles sequence.
2. Set the configurations in the `config.json` file, an example is provided below:
```
{
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
        "MolWt", "MolLogP", "TPSA"
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
```

## Usage
Install the package
```
pip install mol_eval
```
Run the evaluation