from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import Descriptors

from mol_eval.data_loader import DataLoader
from mol_eval.enums import ValidationLabel
from mol_eval.evaluator import MolEvaluator
from mol_eval.schemas import ConfigSchema


@pytest.mark.parametrize(
    "fake_smiles, original_smiles, expected_filtered, test_case_description",
    [
        (
            # Partial overlap case
            pd.DataFrame({"smiles": ["C=O", "C1=CC=CC=C1", "C1CCNCC1", "CCO"]}),
            pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]}),
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            "Partial overlap: Remove overlapping SMILES",
        ),
        (
            # No overlap case
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]}),
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            "No overlap: Keep all SMILES",
        ),
        (
            # Complete overlap case
            pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]}),
            pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]}),
            pd.DataFrame({"smiles": []}),
            "Complete overlap: All SMILES removed",
        ),
        (
            # Empty fake SMILES case
            pd.DataFrame({"smiles": []}),
            pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]}),
            pd.DataFrame({"smiles": []}),
            "Empty fake SMILES: Return empty DataFrame",
        ),
        (
            # Empty original SMILES case
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            pd.DataFrame({"smiles": []}),
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            "Empty original SMILES: Keep all fake SMILES",
        ),
    ],
)
def test_remove_existing_valid_data(
    fake_smiles, original_smiles, expected_filtered, test_case_description
):
    """
    Test the remove_existing method with different valid cases.

    Args:
        fake_smiles (pd.DataFrame): Fake SMILES data.
        original_smiles (pd.DataFrame): Original SMILES data.
        expected_filtered (pd.DataFrame): Expected result after filtering.
        test_case_description (str): Description of the test case.
    """
    evaluator = MolEvaluator()

    # Ensure both DataFrames have the same dtype and column name
    fake_smiles["smiles"] = fake_smiles["smiles"].astype(str)
    original_smiles["smiles"] = original_smiles["smiles"].astype(str)

    # Call the method to filter the DataFrame
    result = evaluator.remove_existing(
        fake_smiles_df=fake_smiles, original_smiles_df=original_smiles
    )

    # Ensure 'smiles' column is treated as string
    expected_filtered["smiles"] = expected_filtered["smiles"].astype(str)

    # Reset index for both result and expected_filtered to avoid index mismatch
    result_reset = result.reset_index(drop=True)
    expected_filtered_reset = expected_filtered.reset_index(drop=True)

    # Assert column names are identical
    assert (
        result_reset.columns.tolist() == expected_filtered_reset.columns.tolist()
    ), f"Failed in: {test_case_description}"

    # Perform DataFrame comparison
    pd.testing.assert_frame_equal(
        result_reset, expected_filtered_reset, check_dtype=True
    )

    # Additional assertion for clarity
    assert isinstance(result_reset, pd.DataFrame), f"Failed in: {test_case_description}"


@pytest.mark.parametrize(
    "fake_smile, real_smiles, threshold, expected_similarity, expected_max_similarity, expected_similar_sequences, "
    "test_case_description",
    [
        (
            "CCO",  # Identical SMILES
            ["CCO", "CCO", "CCO"],  # Real SMILES
            0.5,  # Threshold
            True,  # Expected similarity
            1.0,  # Expected max similarity
            ["CCO", "CCO", "CCO"],  # Expected similar sequences
            "Identical SMILES: All should match with similarity 1.0",
        ),
        (
            "CCO",  # Fake SMILES
            ["CCO", "C1=CC=CC=C1", "CCO"],  # Real SMILES
            0.5,  # Threshold
            True,  # Expected similarity
            1.0,  # Expected max similarity
            ["CCO", "CCO"],  # Expected similar sequences
            "Partial overlap: Some real SMILES should match",
        ),
        (
            "C=O",  # Fake SMILES
            ["CCO", "C1=CC=CC=C1"],  # Real SMILES
            0.8,  # Threshold
            False,  # Expected similarity
            0.0,  # Expected max similarity
            [],  # Expected similar sequences
            "No overlap: No real SMILES should match",
        ),
        (
            "CCO",  # Fake SMILES
            [],  # Empty list of real SMILES
            0.5,  # Threshold
            False,  # Expected similarity
            0.0,  # Expected max similarity
            [],  # Expected similar sequences
            "Empty real SMILES: No matches should occur",
        ),
        (
            "C1CCNCC1",  # Fake SMILES
            ["C1CCNCC1", "CCO", "C1CCNCC1"],  # Real SMILES
            0.5,  # Threshold
            True,  # Expected similarity
            1.0,  # Expected max similarity
            ["C1CCNCC1", "C1CCNCC1"],  # Expected similar sequences
            "Exact match: Fake SMILES matches multiple real SMILES",
        ),
        (
            "C1CCNCC1",  # Fake SMILES
            ["C1CCNCC1", "C1=CC=CC=C1"],  # Real SMILES
            0.7,  # Threshold
            True,  # Expected similarity
            1.0,  # Expected max similarity
            ["C1CCNCC1"],  # Expected similar sequences
            "Threshold match: Only high similarity matches",
        ),
        (
            "C1CCNCC1",  # Fake SMILES
            ["C1=CC=CC=C1"],  # Real SMILES
            0.8,  # Threshold
            False,  # Expected similarity
            0.0,  # Expected max similarity
            [],  # Expected similar sequences
            "No similarity above threshold: No matches",
        ),
    ],
)
def test_compute_similarity(
    fake_smile,
    real_smiles,
    threshold,
    expected_similarity,
    expected_max_similarity,
    expected_similar_sequences,
    test_case_description,
):
    """
    Test the _compute_similarity method with different valid cases using parameterized inputs.

    Args:
        fake_smile (str): Fake SMILES data.
        real_smiles (list): List of real SMILES data.
        threshold (float): Similarity threshold for considering a match.
        expected_similarity (bool): Expected similarity outcome.
        expected_max_similarity (float): Expected maximum similarity.
        expected_similar_sequences (list): Expected list of similar sequences.
        test_case_description (str): Description of the test case.
    """
    evaluator = MolEvaluator()

    # Call the method to compute similarity
    result = evaluator._compute_similarity(fake_smile, real_smiles, threshold)

    # Check the 'similar' value
    assert (
        result["similar"] == expected_similarity
    ), f"Failed in: {test_case_description}"

    # Check the 'max_similarity' value
    assert (
        result["max_similarity"] == expected_max_similarity
    ), f"Failed in: {test_case_description}"

    # Use np.array_equal for array comparison to avoid ambiguity
    assert np.array_equal(
        result["most_similar_sequences"], expected_similar_sequences
    ), f"Failed in: {test_case_description}"

    # Additional assertions for clarity
    assert isinstance(result, dict), f"Failed in: {test_case_description}"
    assert (
        "similar" in result
        and "max_similarity" in result
        and "most_similar_sequences" in result
    ), f"Failed in: {test_case_description}"


@pytest.mark.parametrize(
    "fake_smiles_df, original_smiles_df, threshold, expected_filtered_df, test_case_description",
    [
        (
            # Scenario: Some fake SMILES match with real SMILES
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            pd.DataFrame(
                {
                    "smiles": ["C1=CC=CC=C1", "CCO"],
                    "cmpd_name": ["compound1", "compound2"],
                }
            ),
            0.5,
            pd.DataFrame(
                {
                    "smiles": ["C=O", "C1CCNCC1"],
                    "similar": [False, True],
                    "max_similarity": [0.0, 0.9],
                    "most_similar_sequences": [[], ["C1=CC=CC=C1"]],
                    "matching_cmpd_names": [[], ["compound1"]],
                }
            ),
            "Scenario: Some fake SMILES match with real SMILES",
        ),
        (
            # Scenario: No fake SMILES match with real SMILES
            pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
            pd.DataFrame(
                {
                    "smiles": ["CCO", "C1=CC=CC=C1"],
                    "cmpd_name": ["compound1", "compound2"],
                }
            ),
            0.5,
            pd.DataFrame(
                {
                    "smiles": ["C=O", "C1CCNCC1"],
                    "similar": [False, False],
                    "max_similarity": [0.0, 0.0],
                    "most_similar_sequences": [[], []],
                    "matching_cmpd_names": [[], []],
                }
            ),
            "Scenario: No fake SMILES match with real SMILES",
        ),
        (
            # Scenario: Empty fake SMILES
            pd.DataFrame({"smiles": []}),
            pd.DataFrame({"smiles": ["C1=CC=CC=C1"], "cmpd_name": ["compound1"]}),
            0.5,
            pd.DataFrame(
                {
                    "smiles": [],
                    "similar": [],
                    "max_similarity": [],
                    "most_similar_sequences": [],
                    "matching_cmpd_names": [],
                }
            ),
            "Scenario: Empty fake SMILES, return empty DataFrame",
        ),
    ],
)
def test_add_levenshtein_similarity(
    fake_smiles_df,
    original_smiles_df,
    threshold,
    expected_filtered_df,
    test_case_description,
):
    """
    Test the add_levenshtein_similarity method.

    Args:
        fake_smiles_df (pd.DataFrame): Fake SMILES data.
        original_smiles_df (pd.DataFrame): Original SMILES data.
        threshold (float): Similarity threshold for considering a match.
        expected_filtered_df (pd.DataFrame): Expected result after filtering.
        test_case_description (str): Description of the test case.
    """
    evaluator = MolEvaluator()

    # Ensure that 'smiles' columns are strings
    fake_smiles_df["smiles"] = fake_smiles_df["smiles"].astype(str)
    original_smiles_df["smiles"] = original_smiles_df["smiles"].astype(str)

    # Call the method
    result = evaluator.add_levenshtein_similarity(
        fake_smiles_df=fake_smiles_df,
        original_smiles_df=original_smiles_df,
        threshold=threshold,
    )

    # Reset index to avoid index mismatch
    result_reset = result.reset_index(drop=True)
    expected_filtered_df_reset = expected_filtered_df.reset_index(drop=True)

    # Ensure the DataFrames have the same columns
    assert (
        result_reset.columns.tolist() == expected_filtered_df_reset.columns.tolist()
    ), f"Failed in: {test_case_description}"

    # Ensure that the result is a DataFrame
    assert isinstance(result_reset, pd.DataFrame), f"Failed in: {test_case_description}"


@pytest.mark.parametrize(
    "fake_smiles, thresholds, mock_labels, expected_df",
    [
        # Case 1: Standard solubility labels
        (
            pd.DataFrame({"smiles": ["C=O", "C1=CC=CC=C1", "CCO"]}),
            {"VERY_HIGH": -0.5, "HIGH": 0.5, "MODERATE": 1.5, "LOW": 2.5},
            ["VERY_HIGH", "HIGH", "LOW"],  # Mocked return values
            pd.DataFrame(
                {
                    "smiles": ["C=O", "C1=CC=CC=C1", "CCO"],
                    "solubility_label": ["VERY_HIGH", "HIGH", "LOW"],
                }
            ),
        ),
        # Case 2: Handles invalid SMILES
        (
            pd.DataFrame({"smiles": ["C", "INVALID_SMILES", "O"]}),
            {"VERY_HIGH": -0.5, "HIGH": 0.5, "MODERATE": 1.5, "LOW": 2.5},
            ["HIGH", "INVALID", "MODERATE"],  # Mocked return values
            pd.DataFrame(
                {
                    "smiles": ["C", "INVALID_SMILES", "O"],
                    "solubility_label": ["HIGH", "INVALID", "MODERATE"],
                }
            ),
        ),
    ],
)
@patch("mol_eval.evaluator.MolEvaluator._compute_water_solubility_label")
def test_add_solubility_labels(
    mock_compute_label, fake_smiles, thresholds, mock_labels, expected_df
):
    """
    Test the add_solubility_labels method by mocking _compute_water_solubility_label.

    Args:
        mock_compute_label (Mock): Mocked method.
        fake_smiles (pd.DataFrame): Fake SMILES DataFrame.
        thresholds (dict): Dictionary of solubility thresholds.
        mock_labels (list[str]): Mocked return values.
        expected_df (pd.DataFrame): Expected DataFrame with solubility labels.
    """

    # Mock the return values for _compute_water_solubility_label
    mock_compute_label.side_effect = mock_labels

    # Create MolEvaluator instance
    evaluator = MolEvaluator()

    # Call the method
    result_df = evaluator.add_solubility_labels(fake_smiles, thresholds)

    # Reset index before assertion
    pd.testing.assert_frame_equal(
        result_df.reset_index(drop=True), expected_df.reset_index(drop=True)
    )


@pytest.mark.parametrize(
    "fake_smiles, valid_labels, expected_df",
    [
        (
            pd.DataFrame(
                {
                    "smiles": ["C=O", "C1=CC=CC=C1", "CCO"],
                    "solubility_label": ["low", "high", "low"],
                }
            ),
            ["low"],
            pd.DataFrame(
                {
                    "smiles": ["C=O", "CCO"],
                    "solubility_label": ["low", "low"],
                }
            ),
        ),
        (
            pd.DataFrame(
                {
                    "smiles": ["C", "O", "N"],
                    "solubility_label": ["high", "low", "medium"],
                }
            ),
            ["medium", "high"],
            pd.DataFrame(
                {
                    "smiles": ["C", "N"],
                    "solubility_label": ["high", "medium"],
                }
            ),
        ),
    ],
)
def test_filter_by_solubility(fake_smiles, valid_labels, expected_df):
    """
    Test the filter_by_solubility method.

    Args:
        fake_smiles (pd.DataFrame): Input fake SMILES DataFrame.
        valid_labels (list[str]): List of valid solubility labels.
        expected_df (pd.DataFrame): Expected filtered DataFrame.
    """

    # Call the static method under test
    result_df = MolEvaluator.filter_by_solubility(fake_smiles, valid_labels)

    # Reset indices for comparison
    result_df = result_df.reset_index(drop=True)
    expected_df = expected_df.reset_index(drop=True)

    # Assert DataFrame equality
    pd.testing.assert_frame_equal(result_df, expected_df)


@pytest.mark.parametrize(
    "fake_smiles_df, real_smiles_list, mocked_matches, expected_df, description",
    [
        # Test Case 1: Basic substructure matches
        (
            pd.DataFrame({"smiles": ["C1CCCCC1", "CCO", "C1=CC=CC=C1"]}),
            [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("C1=CC=CC=C1")],
            [[1], [1, 2], [2]],  # Mocked substructure matches
            pd.DataFrame(
                {
                    "smiles": ["C1CCCCC1", "CCO", "C1=CC=CC=C1"],
                    "substructure_matches": [[1], [1, 2], [2]],  # Expected matches
                }
            ),
            "Basic case: Add substructure matches",
        ),
    ],
)
@patch.object(MolEvaluator, "_compute_substructure_matches")
def test_compute_substructure_matches(
    mock_compute_matches,
    fake_smiles_df,
    real_smiles_list,
    mocked_matches,
    expected_df,
    description,
):
    """
    Test `compute_substructure_matches` method.
    """
    evaluator = MolEvaluator()

    # Mock `_compute_substructure_matches` to return predefined results
    mock_compute_matches.side_effect = lambda fake_mol, real_mols: mocked_matches.pop(0)

    # Convert real smiles into RDKit Mol objects
    real_mols = [mol for mol in real_smiles_list if mol is not None]

    # Invoke the method under test
    result_df = evaluator.compute_substructure_matches(
        fake_smiles_df,
        pd.DataFrame({"smiles": [Chem.MolToSmiles(mol) for mol in real_mols]}),
    )

    # Validate the output DataFrame
    try:
        pd.testing.assert_frame_equal(result_df, expected_df, check_dtype=True)
        print(f"Test passed: {description}")
    except AssertionError as e:
        pytest.fail(f"Test failed: {description}\n{e}")


@pytest.mark.parametrize(
    "fake_smiles_with_matches, max_matches, expected_filtered, test_case_description",
    [
        (
            # Basic filtering case
            pd.DataFrame(
                {
                    "smiles": ["C1CCCCC1", "CCO", "C1=CC=CC=C1"],
                    "substructure_matches": [[1, 2], [1], []],
                }
            ),
            1,
            pd.DataFrame(
                {
                    "smiles": ["CCO", "C1=CC=CC=C1"],
                    "substructure_matches": [[1], []],
                }
            ),
            "Filter by max matches",
        ),
    ],
)
def test_filter_by_substructure_matches_number(
    fake_smiles_with_matches, max_matches, expected_filtered, test_case_description
):
    """
    Test filter_by_substructure_matches_number.
    """
    result = MolEvaluator.filter_by_substructure_matches_number(
        fake_smiles_df=fake_smiles_with_matches,
        max_substructure_matches=max_matches,
    )

    # Validate the result
    pd.testing.assert_frame_equal(
        result.reset_index(drop=True), expected_filtered.reset_index(drop=True)
    )


@pytest.mark.parametrize(
    "fake_smiles_data, real_smiles_data, thresholds, expected_similarity, expected_score",
    [
        (
            {"smiles": ["CCO", "C1CCCCC1"]},
            {"smiles": ["CCO", "C1CCCCC1"]},
            {"VERY_HIGH": 0.9, "HIGH": 0.7, "MODERATE": 0.5},
            "HIGH",
            0.8,
        ),
        (
            {"smiles": ["CCO", "C1CCCCO"]},
            {"smiles": ["CCO", "C1CCCCC1"]},
            {"VERY_HIGH": 0.9, "HIGH": 0.7, "MODERATE": 0.5},
            "HIGH",
            0.8,
        ),
        # Add more test cases if needed
    ],
)
@patch.object(MolEvaluator, "compute_tanimoto")
def test_add_tanimoto_similarity_score_and_label(
    mock_compute_tanimoto,
    fake_smiles_data,
    real_smiles_data,
    thresholds,
    expected_similarity,
    expected_score,
):
    # Mock compute_tanimoto to return predictable results
    mock_compute_tanimoto.return_value = {
        "fake_smile": "CCO",
        "max_tanimoto_score": 0.8,
        "tanimoto_similarity": expected_similarity,
        "avg_tanimoto": 0.8,
        "avg_dice": 0.8,
        "most_similar_real_mol": None,
    }

    fake_smiles_df = pd.DataFrame(fake_smiles_data)
    real_smiles_df = pd.DataFrame(real_smiles_data)

    # Call the method
    mol_evaluator = MolEvaluator()
    result_df = mol_evaluator.add_tanimoto_similarity_score_and_label(
        fake_smiles_df, real_smiles_df, thresholds
    )

    # Validate the DataFrame has the added columns and expected results
    assert "max_tanimoto_score" in result_df.columns
    assert "tanimoto_similarity" in result_df.columns
    assert result_df["tanimoto_similarity"].iloc[0] == expected_similarity
    assert result_df["max_tanimoto_score"].iloc[0] == expected_score
    assert result_df["avg_tanimoto"].iloc[0] == 0.8
    assert result_df["avg_dice"].iloc[0] == 0.8


@pytest.mark.parametrize(
    "fake_smiles_data,  expected_result, expected_description",
    [
        (
            {"smiles": ["CCO", "C1CCCCC1"]},
            {"smiles": ["CCO", "C1CCCCC1"]},
            "Filtering all valid mols",
        ),
        (
            {"smiles": ["NON_VALID", "FAKE_MOL_SMILES"]},
            {"smiles": []},
            "Filtering non valid mols",
        ),
    ],
)
def test_remove_non_mols(fake_smiles_data, expected_result, expected_description):
    mol_evaluator = MolEvaluator()

    # Call the method
    result = mol_evaluator.remove_non_molecules(pd.DataFrame(fake_smiles_data))

    # Validate the result
    pd.testing.assert_frame_equal(
        result, pd.DataFrame(expected_result), check_dtype=False
    )


def test_special_characters_in_smiles():
    evaluator = MolEvaluator()

    # Test cases with special characters in SMILES
    smiles_df = pd.DataFrame({"smiles": ["C=O", "C1=CC=CC=C1", "C*O"]})

    # Call the method to filter the SMILES
    result = evaluator.remove_existing(
        fake_smiles_df=smiles_df, original_smiles_df=pd.DataFrame({"smiles": []})
    )

    # Ensure it processes correctly (no special characters issues)
    assert not result.empty
    assert "smiles" in result.columns


def test_similarity_with_zero_threshold():
    evaluator = MolEvaluator()

    fake_smile = "C=O"
    real_smiles = ["CCO", "C1=CC=CC=C1", "C=O"]

    # Expected similarity: only identical SMILES should match
    expected_similarity = {
        "similar": True,
        "max_similarity": 1.0,
        "most_similar_sequences": ["C=O"],
    }

    result = evaluator._compute_similarity(fake_smile, real_smiles, 0.0)

    # Assert that only identical SMILES are considered similar
    assert result == expected_similarity


def test_handle_empty_or_null_values():
    evaluator = MolEvaluator()

    # Create a DataFrame with some null values
    fake_smiles_df = pd.DataFrame({"smiles": ["C=O", "CCO"]})
    original_smiles_df = pd.DataFrame({"smiles": ["C=O", "CCO"]})

    result = evaluator.remove_existing(
        fake_smiles_df=fake_smiles_df, original_smiles_df=original_smiles_df
    )

    # Assert that the None value is properly handled
    assert result.shape[0] == 0  # There should be 0 rows after filtering


def test_empty_dataframe_handling():
    evaluator = MolEvaluator()

    empty_df = pd.DataFrame({"smiles": []})

    # Ensure that the method returns an empty DataFrame when input is empty
    result = evaluator.remove_existing(
        fake_smiles_df=empty_df, original_smiles_df=empty_df
    )
    pd.testing.assert_frame_equal(result, empty_df)


def test_add_qed_score():
    # Create a sample DataFrame with SMILES strings
    data = {"smiles": ["CCO", "CC(=O)O", "invalid_smiles"]}
    df = pd.DataFrame(data)

    # Expected QED scores
    expected_qed = [
        Descriptors.qed(Chem.MolFromSmiles("CCO")),
        Descriptors.qed(Chem.MolFromSmiles("CC(=O)O")),
        0.0,
    ]

    # Call the add_qed_score method
    evaluator = MolEvaluator()
    result_df = evaluator.add_qed_score(df, smiles_column_name="smiles")

    # Assert that the QED scores are as expected
    for i, expected in enumerate(expected_qed):
        if expected is None:
            assert pd.isna(result_df.loc[i, "qed"])
        else:
            assert result_df.loc[i, "qed"] == expected


def test_extract_thresholds(test_config_data):
    config = ConfigSchema(**test_config_data)
    config.LEVENSHTEIN_THRESHOLD = 0.5
    config.TANIMOTO_THRESHOLDS = {"HIGH": 0.8, "MODERATE": 0.5}
    config.SOLUBILITY_THRESHOLDS = {"HIGH": 1.0, "LOW": 0.5}
    config.VALID_SOLUBILITY_LABELS = ["HIGH", "MODERATE"]
    config.VALID_TANIMOTO_LABELS = ["HIGH", "MODERATE"]
    config.MAX_SUBSTRUCTURES_MATCHES = 3
    thresholds = MolEvaluator._extract_thresholds(config)
    assert thresholds["levenshtein"] == 0.5
    assert thresholds["tanimoto"] == {"HIGH": 0.8, "MODERATE": 0.5}
    assert thresholds["solubility"] == {"HIGH": 1.0, "LOW": 0.5}
    assert thresholds["valid_solubility"] == ["HIGH", "MODERATE"]
    assert thresholds["valid_tanimoto"] == ["HIGH", "MODERATE"]
    assert thresholds["max_substructures"] == 3


def test_evaluate_non_molecules():
    df = pd.DataFrame({"smiles": ["CCO", "invalid_smiles", "C1=CC=CC=C1"]})
    df["evaluation"] = ""
    result = MolEvaluator._evaluate_non_molecules(df)
    assert result.loc[1, "evaluation"] == ValidationLabel.NON_MOL.value
    assert result.loc[0, "evaluation"] == ""


def test_evaluate_existing_smiles():
    full_df = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"], "evaluation": ""})
    df = full_df.copy()
    real_smiles_list = ["CCO"]
    dl = DataLoader(real_smiles_list=real_smiles_list, fake_smiles_list=["C1=CC=CC=C1"])
    dl.load_smiles()
    evaluator = MolEvaluator()
    result = evaluator._evaluate_existing_smiles(full_df, df, dl)
    assert len(result) == 1
    assert result.iloc[0]["smiles"] == "C1=CC=CC=C1"
    assert full_df.loc[0, "evaluation"] == ValidationLabel.EXISTING.value


def test_evaluate_duplicates():
    full_df = pd.DataFrame({"smiles": ["CCO", "CCO", "C1=CC=CC=C1"], "evaluation": ""})
    df = full_df.copy()
    evaluator = MolEvaluator()
    result = evaluator._evaluate_duplicates(full_df, df)
    assert len(result) == 2
    assert full_df.loc[1, "evaluation"] == ValidationLabel.DUPLICATE.value


def test_evaluate_levenshtein_similarity():
    df = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]})
    real_smiles_list = ["CCO", "C1CCCCC1"]
    dl = DataLoader(real_smiles_list=real_smiles_list, fake_smiles_list=["C1=CC=CC=C1"])
    dl.load_smiles()
    thresholds = {"levenshtein": 0.8}
    evaluator = MolEvaluator()
    result = evaluator._evaluate_levenshtein_similarity(df, dl, thresholds)
    assert "similar" in result.columns
    assert result.loc[0, "similar"]


def test_evaluate_descriptors():
    df = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"]})
    descriptors = ["MolWt", "MolLogP"]
    evaluator = MolEvaluator()
    result = evaluator._evaluate_descriptors(df, descriptors)
    assert "MolWt" in result.columns
    assert "MolLogP" in result.columns


def test_evaluate_solubility(test_config_data):
    full_df = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"], "evaluation": ""})
    df = full_df.copy()
    thresholds = {
        "solubility": test_config_data["SOLUBILITY_THRESHOLDS"],
        "valid_solubility": ["HIGH"],
    }
    evaluator = MolEvaluator()
    result = evaluator._evaluate_solubility(full_df, df, thresholds)
    assert len(result) == 1
    assert result.iloc[0]["smiles"] == "CCO"
    assert full_df.loc[1, "evaluation"] == ValidationLabel.LOW_SOLUBILITY.value


def test_evaluate_substructure_matches():
    full_df = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"], "evaluation": ""})
    df = full_df.copy()
    dl = DataLoader(real_smiles_list=["CCO"], fake_smiles_list=["C1=CC=CC=C1"])
    dl.load_smiles()
    thresholds = {"max_substructures": 1}
    evaluator = MolEvaluator()
    result = evaluator._evaluate_substructure_matches(full_df, df, dl, thresholds)
    assert len(result) == 1
    assert result.iloc[0]["smiles"] == "C1=CC=CC=C1"
    assert full_df.loc[0, "evaluation"] == ValidationLabel.HIGH_SUBSTRUCTURES.value


def test_evaluate_tanimoto_similarity(test_config_data):
    full_df = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"], "evaluation": ""})
    df = full_df.copy()
    real_smiles_list = ["CCO"]
    dl = DataLoader(real_smiles_list=real_smiles_list, fake_smiles_list=["C1=CC=CC=C1"])
    dl.load_smiles()
    thresholds = {
        'tanimoto': test_config_data['TANIMOTO_THRESHOLDS'],
        'valid_tanimoto': ['LOW', 'HIGH'],
    }
    evaluator = MolEvaluator()
    result = evaluator._evaluate_tanimoto_similarity(full_df, df, dl, thresholds)
    assert len(result) == 1
    assert result.iloc[0]["smiles"] != real_smiles_list[0]
