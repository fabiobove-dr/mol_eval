import re

import pandas as pd
import pytest

from mol_eval.data_loader import DataLoader


@pytest.fixture
def sample_csv(tmp_path):
    """Creates a sample CSV file for testing."""
    real_smiles_path = tmp_path / "real_smiles.csv"
    fake_smiles_path = tmp_path / "fake_smiles.csv"

    real_smiles_data = pd.DataFrame({"smiles": ["CCO", "C1=CC=CC=C1", "C1CCCCC1"]})
    fake_smiles_data = pd.DataFrame({"smiles": ["C=O", "C1=CC=CN=C1", "C1CCNCC1"]})

    real_smiles_data.to_csv(real_smiles_path, index=False)
    fake_smiles_data.to_csv(fake_smiles_path, index=False)

    return real_smiles_path, fake_smiles_path


def test_init_from_files(sample_csv):
    """Test that the DataLoader initializes correctly."""
    real_smiles_path, fake_smiles_path = sample_csv
    loader = DataLoader(
        real_smiles_path=str(real_smiles_path), fake_smiles_path=str(fake_smiles_path)
    )

    assert loader.real_smiles_path == str(real_smiles_path)
    assert loader.fake_smiles_path == str(fake_smiles_path)
    assert loader.real_smiles_df is None
    assert loader.fake_smiles_df is None


def test_init_from_list():
    """Test that the DataLoader initializes correctly with lists."""
    real_smiles_list = ["CCO", "C1=CC=CC=C1", "C1CCCCC1"]
    fake_smiles_list = ["C=O", "C1=CC=CN=C1", "C1CCNCC1"]

    loader = DataLoader(
        real_smiles_list=real_smiles_list, fake_smiles_list=fake_smiles_list
    )

    assert loader.real_smiles_list == real_smiles_list
    assert loader.fake_smiles_list == fake_smiles_list
    assert loader.real_smiles_df is None
    assert loader.fake_smiles_df is None


def test_init_from_path_and_list(sample_csv):
    """Test that the DataLoader initializes correctly with both paths and lists."""
    real_smiles_path, fake_smiles_path = sample_csv
    real_smiles_list = ["CCO", "C1=CC=CC=C1", "C1CCCCC1"]
    fake_smiles_list = ["C=O", "C1=CC=CN=C1", "C1CCNCC1"]

    loader = DataLoader(
        real_smiles_path=str(real_smiles_path),
        fake_smiles_path=str(fake_smiles_path),
        real_smiles_list=real_smiles_list,
        fake_smiles_list=fake_smiles_list,
    )

    assert loader.real_smiles_path == str(real_smiles_path)
    assert loader.fake_smiles_path == str(fake_smiles_path)
    assert loader.real_smiles_list == real_smiles_list
    assert loader.fake_smiles_list == fake_smiles_list
    assert loader.real_smiles_df is None
    assert loader.fake_smiles_df is None


@pytest.mark.parametrize(
    "path, should_exist, expected_exception, expected_message",
    [
        (None, False, ValueError, "The path cannot be None."),
        ("non_existent.csv", False, FileNotFoundError, "The file does not exist: "),
    ],
)
def test_validate_path(
    sample_csv, tmp_path, path, should_exist, expected_exception, expected_message
):
    """Test path validation with different scenarios."""
    real_smiles_path, _ = sample_csv
    loader = DataLoader(real_smiles_path=str(real_smiles_path), fake_smiles_path="")

    if not should_exist and path != "non_existent.csv":
        path = None
    elif not should_exist:
        path = str(tmp_path / path)
        expected_message += re.escape(path)

    if path is None or not should_exist:
        with pytest.raises(expected_exception, match=expected_message):
            loader._validate_path(path)
    else:
        loader._validate_path(path)  # Should not raise an exception


@pytest.mark.parametrize(
    "csv_path, expected_rows, expected_columns",
    [
        ("real_smiles.csv", 3, ["smiles"]),
        ("fake_smiles.csv", 3, ["smiles"]),
    ],
)
def test_load_csv(sample_csv, csv_path, expected_rows, expected_columns):
    """Test loading CSV with different scenarios."""
    real_smiles_path, fake_smiles_path = sample_csv
    loader = DataLoader(
        real_smiles_path=str(real_smiles_path), fake_smiles_path=str(fake_smiles_path)
    )

    if csv_path == "real_smiles.csv":
        csv_path = real_smiles_path
    else:
        csv_path = fake_smiles_path

    df = loader.load_csv(str(csv_path))

    assert isinstance(df, pd.DataFrame)
    assert len(df) == expected_rows
    assert all(col in df.columns for col in expected_columns)


@pytest.mark.parametrize(
    "getter, expected_first_smile",
    [
        ("get_real_smiles", "CCO"),
        ("get_fake_smiles", "C=O"),
    ],
)
def test_get_smiles(sample_csv, getter, expected_first_smile):
    """Test retrieving real and fake SMILES strings."""
    real_smiles_path, fake_smiles_path = sample_csv
    loader = DataLoader(
        real_smiles_path=str(real_smiles_path), fake_smiles_path=str(fake_smiles_path)
    )
    loader.load_smiles()
    smiles = getattr(loader, getter)()

    assert isinstance(smiles, list)
    assert len(smiles) == 3
    assert smiles[0] == expected_first_smile
