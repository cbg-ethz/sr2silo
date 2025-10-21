"""Tests for field mapping between SILO and Loculus."""

from __future__ import annotations

from sr2silo.schema.field_mapping import FieldMapping


def test_field_mapping_init_with_default():
    """Test initialization with default mapping for covid."""
    mapping = FieldMapping("covid")

    assert mapping.organism == "covid"
    assert "sample_id" in mapping.get_silo_fields()
    assert "sampleId" in mapping.get_loculus_fields()


def test_field_mapping_init_with_custom():
    """Test initialization with custom mapping."""
    custom_mapping = {
        "field1": "field1Camel",
        "field2": "field2Camel",
    }
    mapping = FieldMapping("custom_organism", custom_mapping)

    assert mapping.organism == "custom_organism"
    assert mapping.get_silo_fields() == ["field1", "field2"]
    assert mapping.get_loculus_fields() == ["field1Camel", "field2Camel"]


def test_default_mapping_covid():
    """Test default mapping for covid organism."""
    mapping = FieldMapping("covid")

    # Check key mappings from existing code
    assert mapping.silo_to_loculus("sample_id") == "sampleId"
    assert mapping.silo_to_loculus("batch_id") == "batchId"
    assert mapping.silo_to_loculus("location_code") == "locationCode"
    assert mapping.silo_to_loculus("location_name") == "locationName"
    assert mapping.silo_to_loculus("sampling_date") == "samplingDate"
    assert mapping.silo_to_loculus("sr2silo_version") == "sr2siloVersion"


def test_default_mapping_unknown_organism():
    """Test default mapping for unknown organism returns empty."""
    mapping = FieldMapping("unknown_organism")

    assert mapping.get_silo_fields() == []
    assert mapping.get_loculus_fields() == []


def test_silo_to_loculus():
    """Test converting SILO field name to Loculus."""
    mapping = FieldMapping("covid")

    assert mapping.silo_to_loculus("sample_id") == "sampleId"
    assert mapping.silo_to_loculus("nonexistent") is None


def test_loculus_to_silo():
    """Test converting Loculus field name to SILO."""
    mapping = FieldMapping("covid")

    assert mapping.loculus_to_silo("sampleId") == "sample_id"
    assert mapping.loculus_to_silo("nonexistent") is None


def test_add_mapping():
    """Test adding a new field mapping."""
    mapping = FieldMapping("covid")
    initial_count = len(mapping.get_silo_fields())

    mapping.add_mapping("new_field", "newField")

    assert len(mapping.get_silo_fields()) == initial_count + 1
    assert mapping.silo_to_loculus("new_field") == "newField"
    assert mapping.loculus_to_silo("newField") == "new_field"


def test_remove_mapping():
    """Test removing a field mapping."""
    mapping = FieldMapping("covid")
    initial_count = len(mapping.get_silo_fields())

    # Remove an existing mapping
    mapping.remove_mapping("sample_id")

    assert len(mapping.get_silo_fields()) == initial_count - 1
    assert mapping.silo_to_loculus("sample_id") is None


def test_validate_against_schemas_valid():
    """Test validation with valid schemas."""
    mapping = FieldMapping("covid")

    silo_fields = [
        "read_id", "sample_id", "batch_id", "location_code",
        "location_name", "sampling_date", "sr2silo_version"
    ]
    loculus_fields = [
        "date", "sampleId", "batchId", "locationCode",
        "locationName", "samplingDate", "sr2siloVersion"
    ]

    validation = mapping.validate_against_schemas(silo_fields, loculus_fields, strict=False)

    assert validation["missing_silo"] == []
    assert validation["missing_loculus"] == []


def test_validate_against_schemas_missing_silo():
    """Test validation when SILO field in mapping doesn't exist in schema."""
    mapping = FieldMapping("covid")

    # SILO schema missing some fields that are in the mapping
    silo_fields = ["read_id", "sample_id"]  # Missing other fields
    loculus_fields = [
        "date", "sampleId", "batchId", "locationCode",
        "locationName", "samplingDate", "sr2siloVersion"
    ]

    validation = mapping.validate_against_schemas(silo_fields, loculus_fields, strict=False)

    # These fields are in mapping but not in silo_fields
    assert "batch_id" in validation["missing_silo"]
    assert "location_code" in validation["missing_silo"]


def test_validate_against_schemas_missing_loculus():
    """Test validation when Loculus field in mapping doesn't exist in schema."""
    mapping = FieldMapping("covid")

    silo_fields = [
        "read_id", "sample_id", "batch_id", "location_code",
        "location_name", "sampling_date", "sr2silo_version"
    ]
    # Loculus schema missing some fields that are in the mapping
    loculus_fields = ["date", "sampleId"]  # Missing other fields

    validation = mapping.validate_against_schemas(silo_fields, loculus_fields, strict=False)

    # These fields are in mapping but not in loculus_fields
    assert "batchId" in validation["missing_loculus"]
    assert "locationCode" in validation["missing_loculus"]


def test_validate_against_schemas_strict():
    """Test strict validation reports unmapped fields."""
    mapping = FieldMapping("covid")

    silo_fields = [
        "read_id", "sample_id", "batch_id", "location_code",
        "location_name", "sampling_date", "sr2silo_version",
        "extra_silo_field"  # Not in mapping
    ]
    loculus_fields = [
        "date", "sampleId", "batchId", "locationCode",
        "locationName", "samplingDate", "sr2siloVersion",
        "extraLoculusField"  # Not in mapping
    ]

    validation = mapping.validate_against_schemas(silo_fields, loculus_fields, strict=True)

    assert "extra_silo_field" in validation["unmapped_silo"]
    assert "extraLoculusField" in validation["unmapped_loculus"]


def test_is_valid_true():
    """Test is_valid returns True for valid mapping."""
    mapping = FieldMapping("covid")

    silo_fields = [
        "read_id", "sample_id", "batch_id", "location_code",
        "location_name", "sampling_date", "sr2silo_version"
    ]
    loculus_fields = [
        "date", "sampleId", "batchId", "locationCode",
        "locationName", "samplingDate", "sr2siloVersion"
    ]

    assert mapping.is_valid(silo_fields, loculus_fields) is True


def test_is_valid_false():
    """Test is_valid returns False for invalid mapping."""
    mapping = FieldMapping("covid")

    # Missing fields in schemas
    silo_fields = ["read_id", "sample_id"]
    loculus_fields = ["date", "sampleId"]

    assert mapping.is_valid(silo_fields, loculus_fields) is False
