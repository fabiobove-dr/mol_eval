import time
from mol_eval.commons import timeout
import pytest
import json
from mol_eval.commons import load_config_file
from mol_eval.schemas import ConfigSchema


def test_timeout_decorator():
    @timeout(max_timeout=1)
    def long_running_function():
        time.sleep(2)  # Simulate a long-running task
        return "Completed"

    result = long_running_function()
    assert result == []  # Expect the default value when timeout occurs


def test_load_config_file_valid(tmp_path, test_config_data):
    # Create a valid JSON file
    config_path = tmp_path / "config.json"
    with open(config_path, "w") as f:
        json.dump(test_config_data, f)

    # Test loading the valid JSON file
    config = load_config_file(str(config_path))
    assert isinstance(config, ConfigSchema)
    assert config.LEVENSHTEIN_THRESHOLD == 0.5


def test_load_config_file_invalid_extension(tmp_path):
    # Create a file with an invalid extension
    invalid_path = tmp_path / "config.txt"
    with open(invalid_path, "w") as f:
        f.write("{}")

    # Test that it raises a ValueError
    with pytest.raises(ValueError, match="File must be a JSON file."):
        load_config_file(str(invalid_path))


def test_load_config_file_missing_file():
    # Test that it raises a FileNotFoundError for a non-existent file
    with pytest.raises(FileNotFoundError, match="File not found:"):
        load_config_file("non_existent.json")
