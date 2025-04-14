from pydantic import BaseModel, Field, DirectoryPath
from typing import List, Dict


class DatasetSchema(BaseModel):
    smiles_column_name: str = "smiles"
    needed_columns: List[str] = ["smiles", "solubility", "tanimoto"]


class ConfigSchema(BaseModel):
    LEVENSHTEIN_THRESHOLD: float = Field(
        ..., gt=0, description="Threshold for Levenshtein similarity."
    )
    TANIMOTO_THRESHOLDS: Dict[str, float] = Field(
        ...,
        description="Tanimoto similarity thresholds for categories like VERY_HIGH, HIGH, MODERATE, LOW.",
    )
    SOLUBILITY_THRESHOLDS: Dict[str, float] = Field(
        ..., description="Thresholds for solubility classification."
    )
    VALID_SOLUBILITY_LABELS: List[str] = Field(
        ..., description="List of valid solubility labels."
    )
    VALID_TANIMOTO_LABELS: List[str] = Field(
        ..., description="List of valid Tanimoto similarity labels."
    )
    MAX_SUBSTRUCTURES_MATCHES: int = Field(
        ..., ge=0, description="Maximum number of substructure matches allowed."
    )
    REPORT_FOLDER: DirectoryPath = Field(
        ..., description="Path to the directory where reports will be saved."
    )
    RELEVANT_DESCRIPTORS: List[str] = Field(
        ..., description="List of relevant descriptors for molecule evaluation."
    )
