
![icon](icon.png)

# Overview | `mol_eval`

`mol_eval` is a tool for evaluating SMILES data, particularly for distinguishing between real and fake SMILES sequences. It uses configurable thresholds and molecular descriptors to assess similarity and other properties such as solubility.

[![Coverage](https://codecov.io/github/fabiobove-dr/mol_eval/coverage.svg?branch=main)](https://codecov.io/gh/fabiobove-dr/mol_eval)
[![PyPI Latest Release](https://img.shields.io/pypi/v/mol_eval.svg)](https://pypi.org/project/mol_eval/)
![Unit Tests](https://github.com/fabiobove-dr/mol_eval/actions/workflows/codecov.yml/badge.svg)
<br>[![Powered by TaccLab](https://img.shields.io/badge/powered%20by-TaccLab-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://tacclab.org)
[![License](https://img.shields.io/github/license/fabiobove-dr/mol_eval.svg)](https://github.com/fabiobove-dr/mol_eval/blob/main/LICENSE)

---
# Features
- **Real vs Fake SMILES Evaluation:** Compare real and synthetic SMILES sequences based on various similarity thresholds.
- **Similarity Metrics:** Uses Levenshtein distance, Tanimoto coefficient, and molecular descriptors for comparison.
- **Configurable Analysis:** Easily tweak similarity thresholds, solubility labels, and molecular descriptors through a configuration file.
- **Reports:** Generate detailed evaluation reports based on the results.

---
# Installation
To install `mol_eval`, you can use `pip`:
```bash
pip install mol_eval
```

---
# Configuration
Before running the tool, you'll need to prepare your dataset and configuration file.
### Step 1: Prepare Your Dataset Files

### `real_data.csv`<br>
This file should contain two columns:
  - `cmpd_name`: The name of the compound.
  - `smile`: The SMILES string representing the molecule.
### `fake_data.csv`<br>
This file should contain one column:
  - `smile`: The SMILES string of synthetic molecules.

### Step 2: Configuration File (config.json)

The configuration file allows you to set various thresholds and other parameters used in the evaluation. Here's an example configuration file:
```json
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
```
### Key Parameters Explained:

- `Thresholds`: Customize similarity and solubility thresholds for better evaluation.
- `Descriptors`: Choose molecular descriptors for evaluation, such as molecular weight (MolWt), logP (MolLogP), and polar surface area (TPSA).
- `Tanimoto` and `Levenshtein`: Fine-tune the thresholds for calculating molecular similarity.
- `Solubility` Labels: Define the solubility categories based on the solubility values.
- `Report Folder`: Define where to save evaluation reports.
---
# Usage
After installing the package and preparing your dataset and configuration file, you can run the evaluation tool via the command line.
Run the Evaluation

Use the following command to evaluate your datasets:
```bash
mol_eval --real_data /path/to/real_data.csv --fake_data /path/to/fake_data.csv --configs /path/to/config.json
```
```bash
usage: mol_eval [-h] --real_data REAL_DATA --fake_data FAKE_DATA --configs CONFIGS

Molecule Evaluator: Evaluate real and fake SMILES data using a configuration file.

options:
  -h, --help            Show this help message and exit.
  --real_data REAL_DATA Path to the real SMILES data file (CSV).
  --fake_data FAKE_DATA Path to the fake SMILES data file (CSV).
  --configs CONFIGS     Path to the configuration JSON file.
```

---
# Report Generation

The tool generates a report in the folder specified by `REPORT_FOLDER` in the configuration file (default is `./report`). 
The report contains detailed information on the evaluation of the SMILES sequences, including:
- similarity metrics
- solubility predictions
- substructure matching
- 2d visualizations of the molecules and their most similar counterparts

---
# Contributing

Contributions are welcome!
Feel free to open issues or submit pull requests. 
Please ensure all tests pass and that the code follows the PEP 8 style guide.
<br>


### (More Information coming soon!)
## pre-commit 
Add pre-commit hook to run tests before commit
```bash
    pre-commit run --all-files
```
---

# License
This project is licensed under the terms of the GNU General Public License, Version 3.
