# Schema Field Mapping

This module provides simple, organism-specific field mapping for metadata submission to Loculus.

## Purpose

The field mapping configuration maps SILO `ReadMetadata` fields (snake_case) to Loculus submission fields (camelCase). This allows the code to support multiple organisms without hardcoding the mapping logic in `submit_to_loculus.py`.

## Usage

```python
from sr2silo.schema import get_field_mapping

# Get field mapping for an organism
mapping = get_field_mapping("covid")

# mapping is a dict:
# {
#     "sample_id": "sampleId",
#     "batch_id": "batchId",
#     "location_code": "locationCode",
#     ...
# }
```

## Adding New Organisms

To add support for a new organism, simply add a new entry to `ORGANISM_FIELD_MAPPINGS` in `field_mapping.py`:

```python
ORGANISM_FIELD_MAPPINGS = {
    "covid": {
        "sample_id": "sampleId",
        "batch_id": "batchId",
        ...
    },
    "mpox": {
        "sample_id": "sampleId",
        "collection_date": "collectionDate",
        ...
    },
}
```

## Design

The mapping is based on the `ReadMetadata` schema in `src/sr2silo/silo_read_schema.py`, which defines the fields present in SILO read data. The Loculus field names come from the Loculus configuration (see [WisePulse config](https://github.com/cbg-ethz/WisePulse/blob/main/ansible/group_vars/all/main.yml)).

This design is intentionally simple:
- No API fetching (requires connectivity anyway)
- No caching (unnecessary complexity)
- Just a config file that defines the mapping
- Easy to extend for new organisms
