from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest


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
def test_remove_existing_valid_data(fake_smiles, original_smiles, expected_filtered, test_case_description):
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
    fake_smiles['smiles'] = fake_smiles['smiles'].astype(str)
    original_smiles['smiles'] = original_smiles['smiles'].astype(str)

    # Call the method to filter the DataFrame
    result = evaluator.remove_existing(fake_smiles_df=fake_smiles, original_smiles_df=original_smiles)

    # Ensure 'smiles' column is treated as string
    expected_filtered['smiles'] = expected_filtered['smiles'].astype(str)

    # Reset index for both result and expected_filtered to avoid index mismatch
    result_reset = result.reset_index(drop=True)
    expected_filtered_reset = expected_filtered.reset_index(drop=True)

    # Assert column names are identical
    assert result_reset.columns.tolist() == expected_filtered_reset.columns.tolist(), f"Failed in: {test_case_description}"

    # Perform DataFrame comparison
    pd.testing.assert_frame_equal(result_reset, expected_filtered_reset, check_dtype=True)

    # Additional assertion for clarity
    assert isinstance(result_reset, pd.DataFrame), f"Failed in: {test_case_description}"


@pytest.mark.parametrize(
    "fake_smile, real_smiles, threshold, expected_similarity, expected_max_similarity, expected_similar_sequences, test_case_description",
    [
        (
                "CCO",  # Identical SMILES
                ["CCO", "CCO", "CCO"],  # Real SMILES
                0.5,  # Threshold
                True,  # Expected similarity
                1.0,  # Expected max similarity
                ["CCO", "CCO", "CCO"],  # Expected similar sequences
                "Identical SMILES: All should match with similarity 1.0"
        ),
        (
                "CCO",  # Fake SMILES
                ["CCO", "C1=CC=CC=C1", "CCO"],  # Real SMILES
                0.5,  # Threshold
                True,  # Expected similarity
                1.0,  # Expected max similarity
                ["CCO", "CCO"],  # Expected similar sequences
                "Partial overlap: Some real SMILES should match"
        ),
        (
                "C=O",  # Fake SMILES
                ["CCO", "C1=CC=CC=C1"],  # Real SMILES
                0.8,  # Threshold
                False,  # Expected similarity
                0.0,  # Expected max similarity
                [],  # Expected similar sequences
                "No overlap: No real SMILES should match"
        ),
        (
                "CCO",  # Fake SMILES
                [],  # Empty list of real SMILES
                0.5,  # Threshold
                False,  # Expected similarity
                0.0,  # Expected max similarity
                [],  # Expected similar sequences
                "Empty real SMILES: No matches should occur"
        ),
        (
                "C1CCNCC1",  # Fake SMILES
                ["C1CCNCC1", "CCO", "C1CCNCC1"],  # Real SMILES
                0.5,  # Threshold
                True,  # Expected similarity
                1.0,  # Expected max similarity
                ["C1CCNCC1", "C1CCNCC1"],  # Expected similar sequences
                "Exact match: Fake SMILES matches multiple real SMILES"
        ),
        (
                "C1CCNCC1",  # Fake SMILES
                ["C1CCNCC1", "C1=CC=CC=C1"],  # Real SMILES
                0.7,  # Threshold
                True,  # Expected similarity
                1.0,  # Expected max similarity
                ["C1CCNCC1"],  # Expected similar sequences
                "Threshold match: Only high similarity matches"
        ),
        (
                "C1CCNCC1",  # Fake SMILES
                ["C1=CC=CC=C1"],  # Real SMILES
                0.8,  # Threshold
                False,  # Expected similarity
                0.0,  # Expected max similarity
                [],  # Expected similar sequences
                "No similarity above threshold: No matches"
        ),
    ]
)
def test_compute_similarity(fake_smile, real_smiles, threshold, expected_similarity, expected_max_similarity,
                            expected_similar_sequences, test_case_description):
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
    assert result["similar"] == expected_similarity, f"Failed in: {test_case_description}"

    # Check the 'max_similarity' value
    assert result["max_similarity"] == expected_max_similarity, f"Failed in: {test_case_description}"

    # Use np.array_equal for array comparison to avoid ambiguity
    assert np.array_equal(result["most_similar_sequences"],
                          expected_similar_sequences), f"Failed in: {test_case_description}"

    # Additional assertions for clarity
    assert isinstance(result, dict), f"Failed in: {test_case_description}"
    assert "similar" in result and "max_similarity" in result and "most_similar_sequences" in result, f"Failed in: {test_case_description}"


@pytest.mark.parametrize(
    "fake_smiles_df, original_smiles_df, threshold, expected_filtered_df, test_case_description",
    [
        (
                # Scenario: Some fake SMILES match with real SMILES
                pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
                pd.DataFrame({"smiles": ["C1=CC=CC=C1", "CCO"], "cmpd_name": ["compound1", "compound2"]}),
                0.5,
                pd.DataFrame({
                    "smiles": ["C=O", "C1CCNCC1"],
                    "similar": [False, True],
                    "max_similarity": [0.0, 0.9],
                    "most_similar_sequences": [[], ["C1=CC=CC=C1"]],
                    "matching_cmpd_names": [[], ["compound1"]]
                }),
                "Scenario: Some fake SMILES match with real SMILES",
        ),
        (
                # Scenario: No fake SMILES match with real SMILES
                pd.DataFrame({"smiles": ["C=O", "C1CCNCC1"]}),
                pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1"], "cmpd_name": ["compound1", "compound2"]}),
                0.5,
                pd.DataFrame({
                    "smiles": ["C=O", "C1CCNCC1"],
                    "similar": [False, False],
                    "max_similarity": [0.0, 0.0],
                    "most_similar_sequences": [[], []],
                    "matching_cmpd_names": [[], []]
                }),
                "Scenario: No fake SMILES match with real SMILES",
        ),
        (
                # Scenario: Empty fake SMILES
                pd.DataFrame({"smiles": []}),
                pd.DataFrame({"smiles": ["C1=CC=CC=C1"], "cmpd_name": ["compound1"]}),
                0.5,
                pd.DataFrame({"smiles": [], "similar": [], "max_similarity": [], "most_similar_sequences": [],
                              "matching_cmpd_names": []}),
                "Scenario: Empty fake SMILES, return empty DataFrame",
        ),
    ],
)
def test_add_levenshtein_similarity(fake_smiles_df, original_smiles_df, threshold, expected_filtered_df,
                                    test_case_description):
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
    fake_smiles_df['smiles'] = fake_smiles_df['smiles'].astype(str)
    original_smiles_df['smiles'] = original_smiles_df['smiles'].astype(str)

    # Call the method
    result = evaluator.add_levenshtein_similarity(fake_smiles_df=fake_smiles_df, original_smiles_df=original_smiles_df,
                                                  threshold=threshold)

    # Reset index to avoid index mismatch
    result_reset = result.reset_index(drop=True)
    expected_filtered_df_reset = expected_filtered_df.reset_index(drop=True)

    # Ensure the DataFrames have the same columns
    assert result_reset.columns.tolist() == expected_filtered_df_reset.columns.tolist(), f"Failed in: {test_case_description}"

    # Ensure that the result is a DataFrame
    assert isinstance(result_reset, pd.DataFrame), f"Failed in: {test_case_description}"


@pytest.mark.parametrize(
    "fake_smiles, relevant_descriptors, mock_descriptor_values, expected_descriptors",
    [
        (
                pd.DataFrame({"smiles": ["C=O", "C1=CC=CC=C1", "CCO"]}),
                ["desc1", "desc2"],
                {"desc1": 3, "desc2": 1},
                pd.DataFrame({
                    "smiles": ["C=O", "C1=CC=CC=C1", "CCO"],
                    "desc1": [3, 3, 3],
                    "desc2": [1, 1, 1],
                }),
        ),
        (
                pd.DataFrame({"smiles": ["C", "O", "N"]}),
                ["desc1", "desc2"],
                {"desc1": 1, "desc2": 0},
                pd.DataFrame({
                    "smiles": ["C", "O", "N"],
                    "desc1": [1, 1, 1],
                    "desc2": [0, 0, 0],
                }),
        ),
    ],
)
@patch("src.mol_evaluator.mol_evaluator.MolEvaluator._compute_descriptors")
def test_describe_fake_smiles(
        mock_compute_descriptors, fake_smiles, relevant_descriptors, mock_descriptor_values, expected_descriptors
):
    """
    Test the describe_fake_smiles method by mocking the _compute_descriptors method.

    Args:
        mock_compute_descriptors (Mock): Mocked _compute_descriptors method.
        fake_smiles (pd.DataFrame): Input fake SMILES DataFrame.
        relevant_descriptors (list[str]): List of descriptors to compute.
        mock_descriptor_values (dict): Mocked descriptor values to return.
        expected_descriptors (pd.DataFrame): Expected DataFrame with computed descriptors.
    """

    # Mock the _compute_descriptors method to return mock_descriptor_values
    mock_compute_descriptors.return_value = mock_descriptor_values

    # Create an instance of MolEvaluator
    evaluator = MolEvaluator()

    # Call the method under test
    result_df = evaluator.describe_fake_smiles(fake_smiles, relevant_descriptors)

    # Reset indices for comparison
    result_df = result_df.reset_index(drop=True)
    expected_descriptors = expected_descriptors.reset_index(drop=True)

    # Assert DataFrame equality
    pd.testing.assert_frame_equal(result_df, expected_descriptors)


@pytest.mark.parametrize(
    "fake_smiles, thresholds, mock_labels, expected_df",
    [
        (
                pd.DataFrame({"smiles": ["C=O", "C1=CC=CC=C1", "CCO"]}),
                {"low": 0.5, "high": 2.0},
                ["low", "high", "low"],
                pd.DataFrame({
                    "smiles": ["C=O", "C1=CC=CC=C1", "CCO"],
                    "solubility_label": ["low", "high", "low"],
                }),
        ),
        (
                pd.DataFrame({"smiles": ["C", "O", "N"]}),
                {"low": 0.2, "high": 1.0},
                ["high", "low", "medium"],
                pd.DataFrame({
                    "smiles": ["C", "O", "N"],
                    "solubility_label": ["high", "low", "medium"],
                }),
        ),
    ],
)
@patch("src.mol_evaluator.mol_evaluator.MolEvaluator._compute_water_solubility_label")
def test_add_solubility_labels(
        mock_compute_label, fake_smiles, thresholds, mock_labels, expected_df
):
    """
    Test the add_solubility_labels method by mocking _compute_water_solubility_label.

    Args:
        mock_compute_label (Mock): Mocked _compute_water_solubility_label method.
        fake_smiles (pd.DataFrame): Input fake SMILES DataFrame.
        thresholds (dict): Thresholds for solubility labels.
        mock_labels (list[str]): Mocked labels to return.
        expected_df (pd.DataFrame): Expected DataFrame with solubility labels.
    """
    from src.mol_evaluator.mol_evaluator import MolEvaluator

    # Mock return values for each SMILES string
    mock_compute_label.side_effect = mock_labels

    # Create an instance of MolEvaluator
    evaluator = MolEvaluator()

    # Call the method under test
    result_df = evaluator.add_solubility_labels(fake_smiles, thresholds)

    # Reset indices for comparison
    result_df = result_df.reset_index(drop=True)
    expected_df = expected_df.reset_index(drop=True)

    # Assert DataFrame equality
    pd.testing.assert_frame_equal(result_df, expected_df)


@pytest.mark.parametrize(
    "fake_smiles, valid_labels, expected_df",
    [
        (
                pd.DataFrame({
                    "smiles": ["C=O", "C1=CC=CC=C1", "CCO"],
                    "solubility_label": ["low", "high", "low"],
                }),
                ["low"],
                pd.DataFrame({
                    "smiles": ["C=O", "CCO"],
                    "solubility_label": ["low", "low"],
                }),
        ),
        (
                pd.DataFrame({
                    "smiles": ["C", "O", "N"],
                    "solubility_label": ["high", "low", "medium"],
                }),
                ["medium", "high"],
                pd.DataFrame({
                    "smiles": ["C", "N"],
                    "solubility_label": ["high", "medium"],
                }),
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
    from src.mol_evaluator.mol_evaluator import MolEvaluator

    # Call the static method under test
    result_df = MolEvaluator.filter_by_solubility(fake_smiles, valid_labels)

    # Reset indices for comparison
    result_df = result_df.reset_index(drop=True)
    expected_df = expected_df.reset_index(drop=True)

    # Assert DataFrame equality
    pd.testing.assert_frame_equal(result_df, expected_df)


import pytest
import pandas as pd
from unittest.mock import patch
from rdkit import Chem

# Mock class with the methods to test
from src.mol_evaluator.mol_evaluator import MolEvaluator  # Adjust the path as needed


@pytest.mark.parametrize(
    "fake_smiles_df, real_smiles_list, mocked_matches, expected_df, description",
    [
        # Test Case 1: Basic substructure matches
        (
                pd.DataFrame({"smiles": ["C1CCCCC1", "CCO", "C1=CC=CC=C1"]}),
                [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("C1=CC=CC=C1")],
                [[1], [1, 2], [2]],  # Mocked substructure matches
                pd.DataFrame({
                    "smiles": ["C1CCCCC1", "CCO", "C1=CC=CC=C1"],
                    "substructure_matches": [[1], [1, 2], [2]],  # Expected matches
                }),
                "Basic case: Add substructure matches",
        ),
    ],
)
@patch.object(MolEvaluator, '_compute_substructure_matches')
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
    result_df = evaluator.compute_substructure_matches(fake_smiles_df, pd.DataFrame(
        {"smiles": [Chem.MolToSmiles(mol) for mol in real_mols]}))

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
                pd.DataFrame({
                    "smiles": ["C1CCCCC1", "CCO", "C1=CC=CC=C1"],
                    "substructure_matches": [[1, 2], [1], []],
                }),
                1,
                pd.DataFrame({
                    "smiles": ["CCO", "C1=CC=CC=C1"],
                    "substructure_matches": [[1], []],
                }),
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
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected_filtered.reset_index(drop=True))


@pytest.mark.parametrize(
    "fake_smiles_data, real_smiles_data, thresholds, expected_similarity, expected_score",
    [
        (
                {'smiles': ['CCO', 'C1CCCCC1']},
                {'smiles': ['CCO', 'C1CCCCC1']},
                {'VERY_HIGH': 0.9, 'HIGH': 0.7, 'MODERATE': 0.5},
                'HIGH', 0.8
        ),
        (
                {'smiles': ['CCO', 'C1CCCCO']},
                {'smiles': ['CCO', 'C1CCCCC1']},
                {'VERY_HIGH': 0.9, 'HIGH': 0.7, 'MODERATE': 0.5},
                'HIGH', 0.8
        ),
        # Add more test cases if needed
    ]
)
@patch.object(MolEvaluator, 'compute_tanimoto')
def test_add_tanimoto_similarity_score_and_label(mock_compute_tanimoto, fake_smiles_data, real_smiles_data, thresholds,
                                                 expected_similarity, expected_score):
    # Mock compute_tanimoto to return predictable results
    mock_compute_tanimoto.return_value = {
        'fake_smile': 'CCO',
        'max_tanimoto_score': 0.8,
        'tanimoto_similarity': expected_similarity,
        'avg_tanimoto': 0.8,
        'avg_dice': 0.8,
        'most_similar_real_mol': None
    }

    fake_smiles_df = pd.DataFrame(fake_smiles_data)
    real_smiles_df = pd.DataFrame(real_smiles_data)

    # Call the method
    mol_evaluator = MolEvaluator()
    result_df = mol_evaluator.add_tanimoto_similarity_score_and_label(fake_smiles_df, real_smiles_df, thresholds)

    # Validate the DataFrame has the added columns and expected results
    assert 'max_tanimoto_score' in result_df.columns
    assert 'tanimoto_similarity' in result_df.columns
    assert result_df['tanimoto_similarity'].iloc[0] == expected_similarity
    assert result_df['max_tanimoto_score'].iloc[0] == expected_score
    assert result_df['avg_tanimoto'].iloc[0] == 0.8
    assert result_df['avg_dice'].iloc[0] == 0.8


@pytest.mark.parametrize(
    "fake_smiles_data,  expected_result, expected_description",
    [
        (
                {'smiles': ['CCO', 'C1CCCCC1']},
                {'smiles': ['CCO', 'C1CCCCC1']},
                'Filtering all valid mols',
        ),
        (
                {'smiles': ['NON_VALID', 'FAKE_MOL_SMILES']},
                {'smiles': []},
                'Filtering non valid mols',
        ),
    ]
)
def test_remove_non_mols(fake_smiles_data, expected_result, expected_description):
    mol_evaluator = MolEvaluator()

    # Call the method
    result = mol_evaluator.remove_non_molecules(pd.DataFrame(fake_smiles_data))

    # Validate the result
    pd.testing.assert_frame_equal(result, pd.DataFrame(expected_result), check_dtype=False)


if __name__ == "__main__":
    pytest.main()
