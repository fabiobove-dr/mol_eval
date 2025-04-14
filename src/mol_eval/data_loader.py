import os
from typing import List

import pandas as pd
from mol_eval.enums import SmilesData


class DataLoader:
    """A class for loading and managing SMILES data.

    This class provides methods to load and retrieve real and fake SMILES data
    from CSV files or directly from lists.

    Attributes:
        real_smiles_path (str): Path to the CSV file containing real SMILES data.
        fake_smiles_path (str): Path to the CSV file containing fake SMILES data.
        real_smiles_list (pd.DataFrame): Loaded DataFrame of real SMILES data.
        fake_smiles_df (pd.DataFrame): Loaded DataFrame of fake SMILES data.
    """

    def __init__(
        self,
        real_smiles_path: str = None,
        fake_smiles_path: str = None,
        real_smiles_list: list = None,
        fake_smiles_list: list = None,
        smiles_column_name: str = SmilesData.smiles_column_name.value,
    ):
        """Initializes the DataLoader with paths to SMILES data or lists.

        Args:
            real_smiles_path (str, optional): Path to the real SMILES data CSV file.
            fake_smiles_path (str, optional): Path to the fake SMILES data CSV file.
            real_smiles_list (list, optional): List of real SMILES strings.
            fake_smiles_list (list, optional): List of fake SMILES strings.
            smiles_column_name (str, optional): Column name for SMILES strings in the DataFrame.
        """
        self.smiles_column_name = smiles_column_name
        self.real_smiles_path = real_smiles_path
        self.fake_smiles_path = fake_smiles_path
        self.real_smiles_list = real_smiles_list
        self.fake_smiles_list = fake_smiles_list
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
        """Loads both real and fake SMILES data into DataFrames.

        The data is loaded from lists if provided, otherwise from CSV files.
        """
        if self.real_smiles_list is not None:
            self.real_smiles_df = pd.DataFrame(
                {self.smiles_column_name: self.real_smiles_list}
            )
        elif self.real_smiles_path is not None:
            self.real_smiles_df = self.load_csv(self.real_smiles_path)
        else:
            raise ValueError(
                "Either real_smiles_list or real_smiles_path must be provided."
            )

        if self.fake_smiles_list is not None:
            self.fake_smiles_df = pd.DataFrame(
                {self.smiles_column_name: self.fake_smiles_list}
            )
        elif self.fake_smiles_path is not None:
            self.fake_smiles_df = self.load_csv(self.fake_smiles_path)
        else:
            raise ValueError(
                "Either fake_smiles_list or fake_smiles_path must be provided."
            )

    def get_real_smiles(self) -> List[str]:
        """Retrieves a list of real SMILES strings.

        Returns:
            List[str]: A list of SMILES strings from the real SMILES data.

        Raises:
            ValueError: If the real SMILES data has not been loaded.
        """
        if self.real_smiles_df is None:
            raise ValueError("Real SMILES data has not been loaded yet.")
        return self.real_smiles_df[self.smiles_column_name].tolist()

    def get_fake_smiles(self) -> List[str]:
        """Retrieves a list of fake SMILES strings.

        Returns:
            List[str]: A list of SMILES strings from the fake SMILES data.

        Raises:
            ValueError: If the fake SMILES data has not been loaded.
        """
        if self.fake_smiles_df is None:
            raise ValueError("Fake SMILES data has not been loaded yet.")
        return self.fake_smiles_df[self.smiles_column_name].tolist()
