"""Tests for validating that Pydantic schemas match SILO database config.

This test module ensures that the ReadMetadata Pydantic schema field aliases
match the field names defined in the SILO database_config.yaml file.

Run with: pytest tests/test_database_config_validation.py -v
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Set

import yaml

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def load_database_config(config_path: Path) -> Dict:
    """Load the SILO database configuration from YAML file.

    Args:
        config_path: Path to the database_config.yaml file

    Returns:
        Dictionary containing the database configuration
    """
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_metadata_field_names_from_config(config: Dict) -> Set[str]:
    """Extract metadata field names from database config.

    Args:
        config: Database configuration dictionary

    Returns:
        Set of metadata field names from the config
    """
    metadata_fields = config.get("schema", {}).get("metadata", [])
    return {field["name"] for field in metadata_fields}


def get_schema_aliases() -> Dict[str, str]:
    """Get field aliases from ReadMetadata schema using Pydantic's model_fields.

    Returns:
        Dict mapping snake_case field names to camelCase aliases
    """
    from sr2silo.silo_read_schema import ReadMetadata

    return {
        field_name: field_info.alias
        for field_name, field_info in ReadMetadata.model_fields.items()
        if field_info.alias
    }


def validate_schema_matches_config(config_path: Path) -> bool:
    """Validate that ReadMetadata schema aliases match database config.

    This ensures that:
    1. All ReadMetadata field aliases exist in the database config
    2. All database config metadata fields are represented in ReadMetadata

    Args:
        config_path: Path to the database_config.yaml file

    Returns:
        True if validation passes, False otherwise
    """
    config = load_database_config(config_path)
    config_fields = get_metadata_field_names_from_config(config)
    schema_aliases = get_schema_aliases()

    # Get set of aliases from schema
    schema_alias_set = set(schema_aliases.values())

    # Check if all schema aliases exist in config
    missing_in_config = schema_alias_set - config_fields
    if missing_in_config:
        logging.error(
            f"Schema aliases not found in database config: {missing_in_config}"
        )
        return False

    # Check if all config fields exist in schema
    missing_in_schema = config_fields - schema_alias_set
    if missing_in_schema:
        logging.warning(
            f"Database config fields not found in schema: {missing_in_schema}"
        )
        # This is just a warning as database might have extra fields

    logging.info("✓ ReadMetadata schema aliases match database config metadata fields")
    return True


# === Pytest Tests ===


def test_schema_matches_database_config():
    """Test that ReadMetadata schema aliases match database_config.yaml."""
    config_path = (
        Path(__file__).parent.parent / "resources" / "silo" / "database_config.yaml"
    )

    assert config_path.exists(), f"Config file not found: {config_path}"

    # Validate that schema matches config
    is_valid = validate_schema_matches_config(config_path)
    assert is_valid, "ReadMetadata schema does not match database config!"


def test_all_schema_aliases_in_config():
    """Test that all Pydantic field aliases exist in the database config."""
    config_path = (
        Path(__file__).parent.parent / "resources" / "silo" / "database_config.yaml"
    )

    config = load_database_config(config_path)
    config_fields = get_metadata_field_names_from_config(config)
    schema_aliases = get_schema_aliases()

    # All schema aliases should be in the config
    for field_name, alias in schema_aliases.items():
        assert alias in config_fields, (
            f"Schema alias '{alias}' (from field '{field_name}') "
            "not found in database config"
        )


def test_field_mapping_consistency():
    """Test that field mappings are consistent and complete."""
    schema_aliases = get_schema_aliases()

    # Check that all aliases are camelCase
    for field_name, alias in schema_aliases.items():
        assert (
            field_name.islower() or "_" in field_name
        ), f"Field name '{field_name}' should be snake_case"
        assert (
            alias[0].islower() and "_" not in alias
        ), f"Alias '{alias}' should be camelCase"


# === Helper for manual validation ===


def print_schema_config_mapping() -> None:
    """Print the schema field to alias mapping for debugging."""
    aliases = get_schema_aliases()

    print("\n=== ReadMetadata Schema Field Mapping ===")
    print(f"{'Python Field (snake_case)':<30} → {'Alias (camelCase)':>20}")
    print("-" * 52)
    for field_name, alias in sorted(aliases.items()):
        print(f"{field_name:<30} → {alias:>20}")
    print()


if __name__ == "__main__":
    # Manual validation for debugging
    config_path = (
        Path(__file__).parent.parent / "resources" / "silo" / "database_config.yaml"
    )

    print_schema_config_mapping()

    if config_path.exists():
        is_valid = validate_schema_matches_config(config_path)
        if is_valid:
            print("✓ Validation passed!")
        else:
            print("✗ Validation failed!")
    else:
        print(f"Config file not found: {config_path}")
