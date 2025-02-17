from enum import Enum


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
