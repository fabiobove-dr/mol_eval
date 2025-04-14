from mol_eval.enums import SmilesData


def test_molecules_schema():
    expected_value = "smiles"
    assert SmilesData.smiles_column_name.value == expected_value
