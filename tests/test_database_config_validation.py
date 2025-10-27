"""Tests for database configuration validation."""

from __future__ import annotations

from pathlib import Path

import yaml

from sr2silo.silo_read_schema import ReadMetadata

DATABASE_CONFIG = Path("resources/silo/database_config.yaml")


def test_valid_database_config_file():
    """Validates that the schema of the database
    config file matches the ReadMetadata model
    at least in the naming of the fields
    """
    with open(DATABASE_CONFIG, "r") as f:
        config = yaml.safe_load(f)
    metadata_schema = config.get("schema").get("metadata")
    assert metadata_schema is not None, "Missing metadata section in database config"

    # get the first field of each item as the names of the fields
    # get the second field of each item as the type of the fields
    field_names = [field.get("name") for field in metadata_schema]

    pydantic_schema = ReadMetadata.model_json_schema()
    assert pydantic_schema is not None, "ReadMetadata model schema is None"
    pydantic_names = pydantic_schema.get("properties")
    assert pydantic_names is not None, (
        "Missing properties section in ReadMetadata model schema"
    )

    assert sorted(field_names) == sorted(pydantic_names.keys()), (
        "Field names do not match"
    )
