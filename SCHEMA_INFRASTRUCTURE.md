# Schema Infrastructure Implementation

This document describes the organism-specific schema infrastructure implementation for sr2silo.

## Overview

The schema infrastructure provides a robust, future-proof foundation for managing organism-specific configurations in sr2silo. It handles:

1. **SILO Schema Management**: Fetches and caches database schema from the SILO/LAPIS API
2. **Loculus Schema Management**: Loads and parses Loculus configuration (locally stored)
3. **Field Mapping**: Maps between SILO fields (snake_case) and Loculus fields (camelCase)
4. **Validation**: Ensures completeness and correctness of configuration
5. **Metadata Transformation**: Converts SILO metadata to Loculus submission format

## Architecture

```
src/sr2silo/schema/
├── __init__.py              # Module exports
├── silo_schema.py           # SILO schema fetching with caching
├── loculus_schema.py        # Loculus config loading
├── field_mapping.py         # Organism-specific field mappings
├── organism_config.py       # Unified configuration manager
└── README.md                # API documentation

tests/schema/
├── __init__.py
├── test_silo_schema.py      # 14 tests for SILO schema
├── test_loculus_schema.py   # 11 tests for Loculus schema
├── test_field_mapping.py    # 14 tests for field mapping
└── test_organism_config.py  # 11 tests for organism config

examples/
├── schema_usage.py          # Example with live API (requires network)
└── schema_usage_mock.py     # Example with mocked data (no network)
```

## Core Components

### 1. SiloSchema (`silo_schema.py`)

**Purpose**: Fetch SILO database configuration from API with caching support.

**Key Features**:
- Fetches schema from `/info/database-config` endpoint
- Caches schema locally in `~/.sr2silo/cache/silo_schema_{organism}.json`
- Automatic fallback to cache when API is unavailable
- Extracts metadata fields, nucleotide sequences, and genes

**Usage**:
```python
from sr2silo.schema import SiloSchema

silo = SiloSchema("https://lapis.wasap.genspectrum.org", "covid")
schema = silo.fetch_schema(use_cache=True)
metadata_fields = silo.get_metadata_fields()
```

### 2. LoculusSchema (`loculus_schema.py`)

**Purpose**: Load and manage Loculus configuration from local sources.

**Key Features**:
- Loads from dictionary or YAML file
- Parses metadata field definitions
- Extracts field types and requirements
- Identifies file upload categories

**Usage**:
```python
from sr2silo.schema import LoculusSchema

# From dict
loculus = LoculusSchema.from_dict("covid", config_dict)

# From YAML file
loculus = LoculusSchema.from_file("covid", Path("config.yml"))

fields = loculus.get_metadata_fields()
```

### 3. FieldMapping (`field_mapping.py`)

**Purpose**: Define and validate field mappings between SILO and Loculus.

**Key Features**:
- Default mappings for supported organisms (currently COVID)
- Bidirectional field name conversion
- Validation against actual schemas
- Support for custom organism mappings

**Default COVID Mapping**:
| SILO (snake_case)  | Loculus (camelCase) |
|--------------------|---------------------|
| sample_id          | sampleId            |
| batch_id           | batchId             |
| location_code      | locationCode        |
| location_name      | locationName        |
| sampling_date      | samplingDate        |
| sr2silo_version    | sr2siloVersion      |

**Usage**:
```python
from sr2silo.schema import FieldMapping

mapping = FieldMapping("covid")
loculus_field = mapping.silo_to_loculus("sample_id")  # Returns "sampleId"
```

### 4. OrganismConfig (`organism_config.py`)

**Purpose**: Unified configuration manager combining all components.

**Key Features**:
- Single entry point for organism configuration
- Comprehensive validation
- Metadata transformation for submission
- Configuration summary logging

**Usage**:
```python
from sr2silo.schema import OrganismConfig

config = OrganismConfig.from_api_and_config(
    organism="covid",
    lapis_url="https://lapis.wasap.genspectrum.org",
    loculus_config=config_dict,
    use_cache=True
)

# Validate
validation = config.validate(strict=False)
if validation["valid"]:
    print("Configuration is valid")

# Transform metadata
loculus_metadata = config.get_metadata_for_submission(silo_metadata)
```

## Validation System

The validation system ensures:

1. **Schema Completeness**:
   - SILO schema is loaded
   - Loculus schema is loaded

2. **Field Mapping Validity**:
   - All mapped SILO fields exist in SILO schema
   - All mapped Loculus fields exist in Loculus schema
   - No orphaned mappings

3. **Required Fields Coverage**:
   - Loculus required fields are covered by mapping
   - Warnings for missing required field mappings

4. **Strict Mode** (optional):
   - Reports unmapped fields in either schema
   - Helps identify missing mappings

## Testing

All tests pass (48 total):

```bash
cd /home/runner/work/sr2silo/sr2silo
pytest tests/schema/ -v
```

Test coverage includes:
- ✓ API fetching with mocking
- ✓ Caching and cache fallback
- ✓ Schema parsing
- ✓ Field mapping validation
- ✓ Organism configuration validation
- ✓ Metadata transformation
- ✓ Error handling

## Examples

### Basic Usage

```python
from sr2silo.schema import OrganismConfig

# Setup
config = OrganismConfig.from_api_and_config(
    organism="covid",
    lapis_url="https://lapis.wasap.genspectrum.org",
    loculus_config=LOCULUS_CONFIG,
    use_cache=True
)

# Validate
validation = config.validate()
assert validation["valid"]

# Transform
silo_metadata = {
    "sample_id": "A1_05_2024_10_08",
    "sampling_date": "2024-10-08",
    # ...
}
loculus_metadata = config.get_metadata_for_submission(silo_metadata)
# Result: {"sampleId": "A1_05_2024_10_08", "samplingDate": "2024-10-08", ...}
```

### Run Examples

```bash
# Mock example (no network required)
python examples/schema_usage_mock.py

# Live example (requires network access)
python examples/schema_usage.py
```

## Adding Support for New Organisms

To add a new organism (e.g., "mpox"):

1. **Define the field mapping** in `src/sr2silo/schema/field_mapping.py`:

```python
def _get_default_mapping(organism: str) -> Dict[str, str]:
    if organism in ["covid", "sc2", "sars-cov-2"]:
        return { ... }
    elif organism == "mpox":
        return {
            "sample_id": "sampleId",
            "collection_date": "collectionDate",
            # Add mpox-specific mappings
        }
    # ...
```

2. **Provide Loculus configuration** for the organism

3. **Use the schema infrastructure**:

```python
config = OrganismConfig.from_api_and_config(
    organism="mpox",
    lapis_url="https://lapis.mpox.example.org",
    loculus_config=MPOX_LOCULUS_CONFIG,
    use_cache=True
)
```

4. **Add tests** for the new organism

## Integration Points

The schema infrastructure is designed to be integrated into:

### 1. Metadata Creation (Future)
Replace manual field mapping in metadata creation with schema-based transformation:

```python
# Before
metadata = {
    "sampleId": data["sample_id"],
    "batchId": data["batch_id"],
    # ... manual mapping
}

# After
metadata = organism_config.get_metadata_for_submission(silo_data)
```

### 2. Submission Pipeline (Future)
Use validated configuration for submissions:

```python
# Validate before submission
validation = organism_config.validate()
if not validation["valid"]:
    raise Exception(f"Invalid configuration: {validation['errors']}")

# Transform and submit
loculus_metadata = organism_config.get_metadata_for_submission(silo_metadata)
```

### 3. CLI Commands (Future)
Add commands for schema inspection:

```bash
sr2silo schema validate --organism covid
sr2silo schema show --organism covid
sr2silo schema cache-clear --organism covid
```

## Design Principles

1. **Organism-First**: All configuration is organism-specific
2. **Validation-First**: Always validate before use
3. **Fail-Safe**: Cache fallback for API failures
4. **Future-Proof**: Extensible for new organisms
5. **Testable**: Comprehensive test coverage
6. **Documented**: Clear API and examples

## Future Enhancements

Potential improvements for future iterations:

1. **Auto-discovery**: Automatically detect available organisms from API
2. **Schema Versioning**: Track and manage schema versions
3. **Migration Support**: Handle schema changes over time
4. **CLI Integration**: Add schema management commands
5. **Configuration Presets**: Provide pre-configured organism settings
6. **Schema Validation**: Validate SILO data against schema before submission

## References

- **Issue**: cbg-ethz/sr2silo#336 (Schema Infrastructure)
- **WisePulse Config**: https://github.com/cbg-ethz/WisePulse/blob/main/ansible/group_vars/all/main.yml
- **SILO API**: https://lapis.wasap.genspectrum.org/swagger-ui/
- **Existing Code**: `src/sr2silo/submit_to_loculus.py` (field mapping reference)

## Questions or Issues?

For questions about the schema infrastructure:
1. Review this document and `src/sr2silo/schema/README.md`
2. Check the examples in `examples/`
3. Review the test suite in `tests/schema/`
4. Consult the inline documentation in source files
