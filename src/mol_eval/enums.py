from enum import Enum


class SmilesData(Enum):
    smiles_column_name: str = "smiles"


class MolType(Enum):
    """Molecule type enum."""

    REAL = "real"
    FAKE = "fake"


class MolWaterSolubilityLabel(Enum):
    """Molecule water solubility label enum."""

    VERY_HIGH = "VERY_HIGH"
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"
    VERY_LOW = "VERY_LOW"


class ValidationLabel(Enum):
    """Validation label enum."""

    EMPTY = ""
    EXISTING = "existing "
    NON_MOL = "non-molecule "
    DUPLICATE = "duplicate "
    INVALID = "low-solutibility "
    HIGH_TANIMOTO = "high-tanimoto "
    HIGH_SUBSTRUCTURES = "high-substructure "
    LOW_SOLUBILITY = "low-solubility "
