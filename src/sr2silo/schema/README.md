# Schema Infrastructure

This module provides organism-specific schema management for sr2silo, handling the mapping between SILO database fields and Loculus submission fields.

## Overview

The schema infrastructure consists of four main components:

1. **SiloSchema** (`silo_schema.py`) - Fetches and caches SILO database configuration from the API
2. **LoculusSchema** (`loculus_schema.py`) - Loads and manages Loculus configuration
3. **FieldMapping** (`field_mapping.py`) - Maps fields between SILO and Loculus formats
4. **OrganismConfig** (`organism_config.py`) - Combines all components with validation

## Quick Start

### Basic Usage

```python
from sr2silo.schema import OrganismConfig

# Define Loculus configuration (or load from file)
loculus_config = {
    "schema": {
        "metadata": [
            {"name": "sampleId", "type": "string"},
            {"name": "samplingDate", "type": "date"},
            # ... more fields
        ]
    }
}

# Create organism configuration
config = OrganismConfig.from_api_and_config(
    organism="covid",
    lapis_url="https://lapis.wasap.genspectrum.org",
    loculus_config=loculus_config,
    use_cache=True
)

# Validate configuration
validation = config.validate()
if validation["valid"]:
    print("✓ Configuration is valid")
else:
    print("✗ Errors:", validation["errors"])

# Transform metadata for submission
silo_metadata = {
    "sample_id": "A1_05_2024_10_08",
    "sampling_date": "2024-10-08",
    # ... more fields
}

loculus_metadata = config.get_metadata_for_submission(silo_metadata)
# Returns: {"sampleId": "A1_05_2024_10_08", "samplingDate": "2024-10-08", ...}
```

## Components

### SiloSchema

Fetches SILO database schema from the API with optional caching.

```python
from sr2silo.schema import SiloSchema

silo_schema = SiloSchema(
    lapis_url="https://lapis.wasap.genspectrum.org",
    organism="covid"
)

# Fetch schema (uses cache if available)
schema = silo_schema.fetch_schema(use_cache=True)

# Get metadata fields
metadata_fields = silo_schema.get_metadata_fields()
# Returns: ["read_id", "sample_id", "batch_id", ...]

# Get genes
genes = silo_schema.get_genes()
# Returns: ["S", "ORF1a", "N", ...]
```

### LoculusSchema

Loads and parses Loculus configuration.

```python
from sr2silo.schema import LoculusSchema

# From dictionary
loculus_schema = LoculusSchema.from_dict("covid", config_dict)

# From YAML file
loculus_schema = LoculusSchema.from_file("covid", Path("config.yml"))

# Get metadata fields
fields = loculus_schema.get_metadata_fields()
# Returns: ["sampleId", "batchId", "samplingDate", ...]

# Get field type
field_type = loculus_schema.get_field_type("sampleId")
# Returns: "string"
```

### FieldMapping

Manages field name mapping between SILO (snake_case) and Loculus (camelCase).

```python
from sr2silo.schema import FieldMapping

# Uses default mapping for organism
mapping = FieldMapping("covid")

# Convert field names
loculus_field = mapping.silo_to_loculus("sample_id")
# Returns: "sampleId"

silo_field = mapping.loculus_to_silo("samplingDate")
# Returns: "sampling_date"

# Add custom mapping
mapping.add_mapping("custom_field", "customField")

# Validate against schemas
validation = mapping.validate_against_schemas(
    silo_fields=["sample_id", "batch_id"],
    loculus_fields=["sampleId", "batchId"],
    strict=True
)
```

### OrganismConfig

Unified configuration manager combining all components.

```python
from sr2silo.schema import OrganismConfig

# Create from API and config
config = OrganismConfig.from_api_and_config(
    organism="covid",
    lapis_url="https://lapis.wasap.genspectrum.org",
    loculus_config=config_dict,
    use_cache=True
)

# Validate everything
validation = config.validate(strict=False)
print(validation["valid"])      # True/False
print(validation["errors"])     # List of errors
print(validation["warnings"])   # List of warnings

# Transform metadata
loculus_metadata = config.get_metadata_for_submission(silo_metadata)

# Log summary
config.log_configuration_summary()
```

## Default Field Mapping (COVID)

The default mapping for COVID/SARS-CoV-2 is:

| SILO Field (snake_case) | Loculus Field (camelCase) |
|------------------------|---------------------------|
| sample_id              | sampleId                  |
| batch_id               | batchId                   |
| location_code          | locationCode              |
| location_name          | locationName              |
| sampling_date          | samplingDate              |
| sr2silo_version        | sr2siloVersion            |

## Validation

The validation system ensures:

1. **Schema Completeness**: Both SILO and Loculus schemas are loaded
2. **Field Mapping Validity**: All mapped fields exist in their respective schemas
3. **No Extra Fields**: Fields in mapping exist in actual schemas
4. **Required Fields**: Loculus required fields are covered (with warnings)

Validation modes:
- **Normal** (`strict=False`): Checks critical errors only
- **Strict** (`strict=True`): Also reports unmapped fields as warnings

## Caching

SILO schemas are cached locally to improve performance and enable offline use:

- Default cache location: `~/.sr2silo/cache/`
- Cache file format: `silo_schema_{organism}.json`
- Automatic fallback to cache if API is unavailable

```python
# Clear cache for an organism
silo_schema.clear_cache()
```

## Adding New Organisms

To add support for a new organism:

1. Define the Loculus configuration
2. Add field mapping in `field_mapping.py` `_get_default_mapping()`
3. Create organism configuration:

```python
def _get_default_mapping(organism: str) -> Dict[str, str]:
    if organism in ["covid", "sc2", "sars-cov-2"]:
        return { ... }
    elif organism == "mpox":
        return {
            "sample_id": "sampleId",
            "collection_date": "collectionDate",
            # ... mpox-specific mappings
        }
    # ...
```

## Testing

Comprehensive test suite in `tests/schema/`:

```bash
pytest tests/schema/ -v
```

Tests cover:
- SILO schema fetching and caching
- Loculus schema loading from dict and file
- Field mapping validation
- Organism configuration validation
- Metadata transformation

## Example

See `examples/schema_usage.py` for a complete working example.

```bash
cd /home/runner/work/sr2silo/sr2silo
python examples/schema_usage.py
```

## Future Integration

This schema infrastructure will be integrated into:
1. Metadata creation logic (validating fields before submission)
2. CLI commands (for schema inspection and validation)
3. Submission pipeline (automatic field transformation)

## References

- WisePulse Loculus Config: https://github.com/cbg-ethz/WisePulse/blob/main/ansible/group_vars/all/main.yml
- SILO API Documentation: https://lapis.wasap.genspectrum.org/swagger-ui/
