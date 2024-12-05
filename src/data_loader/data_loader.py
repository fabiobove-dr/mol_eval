import os
from concurrent.futures import ThreadPoolExecutor
from typing import List

import pandas as pd


class DataLoader:
    """A class for loading and managing SMILES data.

    This class provides methods to load and retrieve real and fake SMILES data
    from CSV files in parallel.

    Attributes:
        real_smiles_path (str): Path to the CSV file containing real SMILES data.
        fake_smiles_path (str): Path to the CSV file containing fake SMILES data.
        real_smiles_df (pd.DataFrame): Loaded DataFrame of real SMILES data.
        fake_smiles_df (pd.DataFrame): Loaded DataFrame of fake SMILES data.
    """

    def __init__(self, real_smiles_path: str, fake_smiles_path: str):
        """Initializes the DataLoader with paths to real and fake SMILES data.

        Args:
            real_smiles_path (str): Path to the real SMILES data CSV file.
            fake_smiles_path (str): Path to the fake SMILES data CSV file.
        """
        self.real_smiles_path = real_smiles_path
        self.fake_smiles_path = fake_smiles_path
        self.real_smiles_df = None
        self.fake_smiles_df = None

    @staticmethod
    def _validate_path(path: str) -> None:
        """Validates that the given file path exists and is not None.

        Args:
            path (str): The file path to validate.

        Raises:
            ValueError: If the path is None.
            FileNotFoundError: If the file does not exist.
        """
        if path is None:
            raise ValueError("The path cannot be None.")
        if not os.path.exists(path):
            raise FileNotFoundError(f"The file does not exist: {path}")

    def load_csv(self, path: str) -> pd.DataFrame:
        """Loads a CSV file into a pandas DataFrame after validating the path.

        Args:
            path (str): Path to the CSV file.

        Returns:
            pd.DataFrame: Loaded DataFrame containing the file's data.
        """
        self._validate_path(path)
        return pd.read_csv(path)

    def load_smiles(self) -> None:
        """Loads both real and fake SMILES data into DataFrames in parallel.

        The data is stored in `real_smiles_df` and `fake_smiles_df` attributes.
        """
        paths = [self.real_smiles_path, self.fake_smiles_path]
        attr_names = ["real_smiles_df", "fake_smiles_df"]

        def _load_and_set_attr(path: str, attr_name: str):
            """Loads a CSV file and assigns it to the corresponding attribute."""
            setattr(self, attr_name, self.load_csv(path))

        with ThreadPoolExecutor() as executor:
            executor.map(_load_and_set_attr, paths, attr_names)

    def get_real_smiles(self) -> List[str]:
        """Retrieves a list of real SMILES strings.

        Returns:
            List[str]: A list of SMILES strings from the real SMILES data.

        Raises:
            ValueError: If the real SMILES data has not been loaded.
        """
        if self.real_smiles_df is None:
            raise ValueError("Real SMILES data has not been loaded yet.")
        return self.real_smiles_df['smiles'].tolist()

    def get_fake_smiles(self) -> List[str]:
        """Retrieves a list of fake SMILES strings.

        Returns:
            List[str]: A list of SMILES strings from the fake SMILES data.

        Raises:
            ValueError: If the fake SMILES data has not been loaded.
        """
        if self.fake_smiles_df is None:
            raise ValueError("Fake SMILES data has not been loaded yet.")
        return self.fake_smiles_df['smiles'].tolist()
