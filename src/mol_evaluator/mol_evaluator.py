import base64
import concurrent.futures
import io
import os

import Levenshtein
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit import RDLogger
from rdkit.Chem import Descriptors, Draw

from commons import timeout

RDLogger.DisableLog('rdApp.*')


class MolEvaluator:
    def __init__(self):
        pass

    @staticmethod
    def remove_duplicates(fake_smiles_df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove duplicate smiles from the fake dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with duplicate smiles removed.
        """
        if not isinstance(fake_smiles_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        if 'smiles' not in fake_smiles_df.columns:
            raise ValueError("DataFrame must contain a 'smiles' column.")

        fake_smiles_df = fake_smiles_df.drop_duplicates(subset=['smiles'])

        return fake_smiles_df

    @staticmethod
    def remove_non_molecules(fake_smiles_df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove non-molecules from the fake dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with non-molecules removed.
        """
        if not isinstance(fake_smiles_df, pd.DataFrame):
            raise TypeError("Input must be a pandas DataFrame.")

        if 'smiles' not in fake_smiles_df.columns:
            raise ValueError("DataFrame must contain a 'smiles' column.")

        for smiles in fake_smiles_df['smiles']:
            if not Chem.MolFromSmiles(smiles):
                fake_smiles_df.drop(fake_smiles_df[fake_smiles_df['smiles'] == smiles].index, inplace=True)

        return fake_smiles_df

    @staticmethod
    def remove_existing(fake_smiles_df: pd.DataFrame, original_smiles_df: pd.DataFrame) -> pd.DataFrame:
        """
        Remove existing smiles from the fake dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            original_smiles_df (pd.DataFrame): Original SMILES dataframe.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with existing SMILES removed.
        """

        if not isinstance(fake_smiles_df, pd.DataFrame) or not isinstance(original_smiles_df, pd.DataFrame):
            raise TypeError("Both inputs must be pandas DataFrames.")

        if 'smiles' not in fake_smiles_df.columns or 'smiles' not in original_smiles_df.columns:
            raise ValueError("Both DataFrames must contain a 'smiles' column.")

        return fake_smiles_df[~fake_smiles_df['smiles'].isin(original_smiles_df['smiles'])]

    @staticmethod
    def _compute_similarity(fake_smile: str, real_smiles: list[str], threshold: float = 0.5) -> dict:
        """
        Compute the similarity between a fake smile and a list of real smiles using Levenshtein distance.

        Args:
            fake_smile (str): Fake SMILES string.
            real_smiles (list): List of real SMILES strings.
            threshold (float): Similarity threshold for considering a match.

        Returns:
            dict: A dictionary containing:
                - `similar`: Boolean indicating if there is a similar SMILES.
                - `max_similarity`: Maximum similarity score.
                - `most_similar_sequences`: List of most similar real SMILES sequences.
        """
        # precompute the similarity scores using Levenshtein distance in a vectorized way
        fake_len = len(fake_smile)
        real_lengths = np.array([len(real_smile) for real_smile in real_smiles])

        # vectorized Levenshtein distance calculation
        similarities = np.array([
            1 - Levenshtein.distance(fake_smile, real_smile) / max(fake_len, real_len)
            for real_smile, real_len in zip(real_smiles, real_lengths)
        ])

        # filter out sequences with similarity below the threshold
        valid_similarities = similarities >= threshold
        similar_sequences = np.array(real_smiles)[valid_similarities]
        valid_similarities = similarities[valid_similarities]

        # if no similar sequences meet the threshold
        similar = valid_similarities.size > 0

        # if there are similar sequences, find the ones with maximum similarity
        if similar:
            max_similarity = np.max(valid_similarities)
            most_similar_sequences = similar_sequences[valid_similarities == max_similarity]
        else:
            max_similarity = 0.0
            most_similar_sequences = []

        return {
            "similar": similar,
            "max_similarity": max_similarity,
            "most_similar_sequences": most_similar_sequences
        }

    def add_levenshtein_similarity(self, fake_smiles_df: pd.DataFrame, original_smiles_df: pd.DataFrame,
                                   threshold: float = 0.5) -> pd.DataFrame:
        """
        Add Levenshtein similarity from the fake dataframe based on similarity to real smiles.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            original_smiles_df (pd.DataFrame): Original SMILES dataframe.
            threshold (float): Similarity threshold for considering a match.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with existing SMILES removed.
        """

        real_smiles_list: list[str] = original_smiles_df['smiles'].tolist()
        fake_smiles_list: list[str] = fake_smiles_df['smiles'].tolist()

        filtered_data = []
        for fake_smiles in fake_smiles_list:
            # Compute similarity using pre-defined method
            similar, max_similarity, most_similar_sequences = self._compute_similarity(
                fake_smiles, real_smiles_list, threshold
            ).values()

            # Ensure that most_similar_sequences is a list, even if it's a single string
            if isinstance(most_similar_sequences, str):
                most_similar_sequences = [most_similar_sequences]

            # Efficiently match cmpd_names using pandas merge
            matching_cmpd_names = self._get_matching_cmpd_names(most_similar_sequences, original_smiles_df)

            # Store the result in a dictionary
            filtered_data.append({
                'smiles': fake_smiles,
                'similar': similar,
                'max_similarity': max_similarity,
                'most_similar_sequences': most_similar_sequences,
                'matching_cmpd_names': matching_cmpd_names
            })

        if not filtered_data:
            filtered_data = [{
                'smiles': "",
                'similar': False,
                'max_similarity': 0.0,
                'most_similar_sequences': [],
                'matching_cmpd_names': []
            }]

        # Create a DataFrame with filtered data
        return pd.DataFrame(filtered_data)

    @staticmethod
    def _get_matching_cmpd_names(most_similar_sequences: list[str], original_smiles_df: pd.DataFrame) -> list[str]:
        """
        Efficiently match compound names for the most similar sequences.

        Args:
            most_similar_sequences (list[str]): List of most similar sequences.
            original_smiles_df (pd.DataFrame): Original SMILES dataframe.

        Returns:
            list: List of matched compound names.
        """
        # Ensure most_similar_sequences is a list, even if it's a single string
        if isinstance(most_similar_sequences, str):
            most_similar_sequences = [most_similar_sequences]

        # Merge the sequences with cmpd_name in the original dataframe
        matched_data = original_smiles_df[original_smiles_df['smiles'].isin(most_similar_sequences)]
        return matched_data['cmpd_name'].tolist()

    @staticmethod
    def _compute_descriptors(fake_smiles: str) -> dict[str, float]:
        """
        Compute the descriptors for a given SMILES string.

        Args:
            fake_smiles (str): SMILES string.

        Returns:
            dict: Dictionary of descriptors and their values.
        """
        mol = Chem.MolFromSmiles(fake_smiles)
        descriptors = {desc[0]: desc[1](mol) for desc in Descriptors.descList}
        return descriptors

    def describe_fake_smiles(self, fake_smiles_df: pd.DataFrame, relevant_descriptors: list[str]) -> pd.DataFrame:
        """
        Add descriptors to the fake SMILES dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            relevant_descriptors (list[str]): List of relevant descriptors.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with descriptors.
        """

        def compute_and_assign_descriptors(row):
            # Compute descriptors for the given SMILES
            descriptors = self._compute_descriptors(row['smiles'])
            # Assign each descriptor to its corresponding column
            for descriptor in relevant_descriptors:
                row[str(descriptor)] = descriptors[descriptor]
            return row

        # Apply the function to compute descriptors row-wise
        return fake_smiles_df.apply(compute_and_assign_descriptors, axis=1)

    @staticmethod
    def _compute_water_solubility_label(smiles: str, thresholds: dict) -> str:
        """
        Compute the water solubility label for a given SMILES string.

        Args:
            smiles (str): SMILES string.

        Returns:
            str: Water solubility label.
        """
        mol = Chem.MolFromSmiles(smiles)
        log_p = Descriptors.MolLogP(mol)
        if log_p < thresholds["VERY_HIGH"]:
            return "VERY_HIGH"
        elif thresholds["VERY_HIGH"] <= log_p < thresholds["HIGH"]:
            return "HIGH"
        elif thresholds["HIGH"] <= log_p < thresholds["MODERATE"]:
            return "MODERATE"
        elif thresholds["MODERATE"] <= log_p < thresholds["LOW"]:
            return "LOW"
        else:
            return "VERY_LOW"

    def add_solubility_labels(self, fake_smiles_df: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
        """
        Add water solubility labels to the fake SMILES dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            thresholds (dict): Dictionary of water solubility thresholds.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with water solubility labels.
        """

        for index, row in fake_smiles_df.iterrows():
            fake_smiles_df.at[index, 'solubility_label'] = self._compute_water_solubility_label(row['smiles'],
                                                                                                thresholds)

        return fake_smiles_df

    @staticmethod
    def filter_by_solubility(fake_smiles_df: pd.DataFrame, valid_labels: list[str]) -> pd.DataFrame:
        """
        Filter the fake SMILES dataframe by water solubility labels.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            valid_labels (list[str]): List of valid water solubility labels.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with filtered water solubility.
        """
        return fake_smiles_df[fake_smiles_df['solubility_label'].isin(valid_labels)]

    @timeout(5.0)
    def _compute_substructure_matches(self, fake_mol: Chem.Mol, real_mols: Chem.Mol) -> list:
        """
        Compute substructure matches for a given SMILES string.
        If the total time exceeds 5 seconds, an exception is raised and an empty list is returned.
        Args:
            fake_mol (Chem.Mol): Fake molecule.
            real_mols (Chem.Mol): Real molecules.
        Returns:
            list: List of matched compound names.
        """
        sub_molecules = []
        for patented_molecule in real_mols:
            try:
                if fake_mol.HasSubstructMatch(patented_molecule):
                    sub_molecules.append(Chem.MolToSmiles(patented_molecule))
                if patented_molecule.HasSubstructMatch(fake_mol):
                    sub_molecules.append(Chem.MolToSmiles(patented_molecule))
            except Exception as _:
                # Handle any issues during the substructure matching process
                sub_molecules.append([])  # Adding empty list in case of error
        return sub_molecules

    def compute_substructure_matches(self, fake_smiles_df: pd.DataFrame, real_smiles_df: pd.DataFrame) -> pd.DataFrame:
        """
        Add substructure matches to the fake SMILES dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            real_smiles_df (pd.DataFrame): Real SMILES dataframe.
        Returns:
            pd.DataFrame: Fake SMILES dataframe with substructure matches.
        """
        # Convert real SMILES to RDKit Molecule objects
        real_mols = [Chem.MolFromSmiles(smiles) for smiles in real_smiles_df['smiles']]

        # Initialize the result column
        fake_smiles_df['substructure_matches'] = None

        for index, fake_smiles in fake_smiles_df.iterrows():
            fake_mol = Chem.MolFromSmiles(fake_smiles['smiles'])
            if not fake_mol:
                fake_smiles_df.at[index, 'substructure_matches'] = []
                continue

            # Compute matches
            try:
                matches = self._compute_substructure_matches(fake_mol, real_mols)
                fake_smiles_df.at[index, 'substructure_matches'] = matches
            except Exception as e:
                fake_smiles_df.at[index, 'substructure_matches'] = []
                print(f"Error computing substructure matches for index {index}: {e}")

        return fake_smiles_df

    @staticmethod
    def filter_by_substructure_matches_number(fake_smiles_df: pd.DataFrame,
                                              max_substructure_matches: int) -> pd.DataFrame:
        """
        Filter the fake SMILES dataframe by substructure matches number.

        Args:
            fake_smiles_df (pd.DataFrame): Fake SMILES dataframe.
            max_substructure_matches (int): Maximum number of substructure matches.

        Returns:
            pd.DataFrame: Fake SMILES dataframe with filtered substructure matches number.
        """
        return fake_smiles_df[fake_smiles_df['substructure_matches'].apply(len) <= max_substructure_matches]

    @staticmethod
    def _visualize_2d(canonical_smiles, as_html=False, is_mol=False):
        """
        Generate a 2D visualization of a molecule from its SMILES string.

        Args:
            canonical_smiles (str): SMILES string of the molecule.
            as_html (bool): If True, return the image as an HTML <img> tag.
            is_mol (bool): If True, interpret `canonical_smiles` as an RDKit molecule.

        Returns:
            Tuple: (Image object or HTML string, RDKit molecule object)
        """
        # Create a molecule object from SMILES string
        mol = Chem.MolFromSmiles(canonical_smiles) if not is_mol else canonical_smiles

        # Check if molecule is valid
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Generate the 2D image of the molecule
        img = Draw.MolToImage(mol, size=(300, 300))

        # Save the image to a bytes buffer
        buf = io.BytesIO()
        img.save(buf, format='PNG')
        buf.seek(0)

        # Encode the image in base64
        img_base64 = base64.b64encode(buf.getvalue()).decode('utf-8')

        # Generate HTML image tag
        img_html = f'<img src="data:image/png;base64,{img_base64}" width="300" height="300">'
        return img_html

    import pandas as pd

    def add_2d_visualizations(self, fake_smiles_df: pd.DataFrame) -> pd.DataFrame:
        """
        Add 2D molecule visualizations to the fake SMILES dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): DataFrame containing a 'smiles' column.

        Returns:
            pd.DataFrame: DataFrame with an added '2d_image' column containing base64-encoded images.
        """
        images = []
        most_sim_img = []
        for index, row in fake_smiles_df.iterrows():
            # Generate the 2D image and convert it to base64
            img = self._visualize_2d(row['smiles'], as_html=False)
            sim_img = self._visualize_2d(row['most_similar_real_mol'], as_html=False)
            images.append(img)
            most_sim_img.append(sim_img)
        fake_smiles_df["2d_image"] = images
        fake_smiles_df["2d_image_most_sim_similes"] = most_sim_img
        return fake_smiles_df

    @staticmethod
    def compute_tanimoto(fake_smile: str, real_fps: dict, thresholds: dict) -> dict:
        """
        Compute Tanimoto and Dice similarity scores for a fake SMILE against real molecules.

        Args:
            fake_smile (str): The fake SMILE string.
            real_fps (dict): Precomputed fingerprints for real molecules.
            thresholds (dict): Tanimoto threshold for similarity score.
        Returns:
            dict: A dictionary with similarity scores, most similar real molecule, and similarity labels.
        """
        fake_mol = Chem.MolFromSmiles(fake_smile)
        if fake_mol is None:
            raise ValueError(f"Invalid SMILE: {fake_smile}")
        fp_fake = Chem.RDKFingerprint(fake_mol)

        sequence_tanimoto_scores, sequence_dice_scores = [], []
        highest_tanimoto_score = 0
        most_similar_real_mol = ''

        for real_smile, fp_real in real_fps.items():
            # Calculate Tanimoto and Dice similarity
            tanimoto_score = DataStructs.TanimotoSimilarity(fp_real, fp_fake)
            dice_score = DataStructs.DiceSimilarity(fp_real, fp_fake)

            sequence_tanimoto_scores.append(tanimoto_score)
            sequence_dice_scores.append(dice_score)

            # Track the most similar real molecule
            if tanimoto_score > highest_tanimoto_score:
                highest_tanimoto_score = tanimoto_score
                most_similar_real_mol = real_smile

        # Calculate average similarity scores
        avg_tanimoto = np.mean(sequence_tanimoto_scores)
        avg_dice = np.mean(sequence_dice_scores)

        # Determine similarity category
        if highest_tanimoto_score >= thresholds['VERY_HIGH']:
            similarity = "VERY HIGH"
        elif highest_tanimoto_score >= thresholds['HIGH']:
            similarity = "HIGH"
        elif highest_tanimoto_score < thresholds['MODERATE']:
            similarity = "LOW"
        else:
            similarity = "MODERATE"

        return {
            'fake_smile': fake_smile,
            'max_tanimoto_score': highest_tanimoto_score,
            'max_dice_score': max(sequence_dice_scores),
            'most_similar_real_mol': most_similar_real_mol,
            'tanimoto_similarity': similarity,
            'avg_tanimoto': avg_tanimoto,
            'avg_dice': avg_dice
        }

    def add_tanimoto_similarity_score_and_label(self, fake_smiles_df: pd.DataFrame,
                                                real_smiles_df: pd.DataFrame, thresholds: dict) -> pd.DataFrame:
        """
        Add Tanimoto similarity scores and labels to the fake SMILES dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): DataFrame containing fake SMILES strings in the 'smiles' column.
            real_smiles_df (pd.DataFrame): DataFrame containing real SMILES strings in the 'smiles' column.
            thresholds (dict): Tanimoto threshold for similarity score.
        Returns:
            pd.DataFrame: Fake SMILES dataframe with added similarity scores and labels.
        """
        # Precompute fingerprints for real molecules
        real_fps = {
            real_smile: Chem.RDKFingerprint(Chem.MolFromSmiles(real_smile))
            for real_smile in real_smiles_df["smiles"]
        }

        results = []
        try:
            with concurrent.futures.ThreadPoolExecutor() as executor:
                for result in executor.map(lambda fake: self.compute_tanimoto(fake, real_fps, thresholds),
                                           fake_smiles_df["smiles"]):
                    results.append(result)
        except Exception as e:
            raise e
        finally:
            executor.shutdown(wait=True)

        # Convert results into a DataFrame and merge
        results_df = pd.DataFrame(results)
        fake_smiles_df = fake_smiles_df.merge(results_df, left_on='smiles', right_on='fake_smile', how='left')
        fake_smiles_df.drop(columns='fake_smile', inplace=True)

        return fake_smiles_df

    @staticmethod
    def filter_by_tanimoto_label(fake_smiles_df: pd.DataFrame, allowed_labels: list) -> pd.DataFrame:
        """
        Filter the fake SMILES dataframe based on allowed Tanimoto similarity labels.

        Args:
            fake_smiles_df (pd.DataFrame): DataFrame containing fake SMILES strings with similarity labels.
            allowed_labels (list): List of allowed similarity labels (e.g., ["LOW", "MODERATE", "HIGH"]).

        Returns:
            pd.DataFrame: Filtered fake SMILES dataframe.
        """
        return fake_smiles_df[fake_smiles_df["tanimoto_similarity"].isin(allowed_labels)]

    @staticmethod
    def create_report(fake_smiles_df: pd.DataFrame, report_folder: str) -> None:
        """
        Create a report for the fake SMILES dataframe.

        Args:
            fake_smiles_df (pd.DataFrame): DataFrame containing fake SMILES strings with similarity labels.
            report_folder (str): Path to the directory where the report will be saved.
        """
        # Create report folder if it doesn't exist
        os.makedirs(report_folder, exist_ok=True)

        # Save report
        fake_smiles_df.to_csv(os.path.join(report_folder, "report.csv"), index=False)
        html_table = fake_smiles_df.to_html(escape=False)  # Set escape=False to allow HTML rendering of images
        with open(os.path.join(report_folder, "report.html"), "w") as f:
            f.write(html_table)
