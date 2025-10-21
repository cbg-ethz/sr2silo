"""Mock example usage of the schema infrastructure (no network required).

This example demonstrates how to use the schema infrastructure with mocked data,
useful for testing and development without network access.
"""

from sr2silo.schema import (
    FieldMapping,
    LoculusSchema,
    OrganismConfig,
    SiloSchema,
)

# Mock SILO schema
MOCK_SILO_SCHEMA = {
    "metadata": [
        {"name": "read_id", "type": "string"},
        {"name": "sample_id", "type": "string"},
        {"name": "batch_id", "type": "string"},
        {"name": "sampling_date", "type": "date"},
        {"name": "location_name", "type": "string"},
        {"name": "location_code", "type": "string"},
        {"name": "sr2silo_version", "type": "string"},
    ],
    "nucleotideSequences": [{"name": "main", "sequence": "ACGT..."}],
    "genes": [
        {"name": "S", "sequence": "MFV..."},
        {"name": "ORF1a", "sequence": "MES..."},
        {"name": "N", "sequence": "MPR..."},
    ],
}

# Mock Loculus configuration (similar to WisePulse)
MOCK_LOCULUS_CONFIG = {
    "schema": {
        "organismName": "SARS-CoV-2",
        "metadata": [
            {"name": "date", "type": "date", "displayName": "Submission Date"},
            {"name": "sampleId", "type": "string", "displayName": "Sample ID"},
            {"name": "batchId", "type": "string", "displayName": "Batch ID"},
            {"name": "locationCode", "type": "string", "displayName": "Location Code"},
            {"name": "locationName", "type": "string", "displayName": "Location"},
            {"name": "samplingDate", "type": "date", "displayName": "Sampling Date"},
            {"name": "sr2siloVersion", "type": "string", "displayName": "sr2silo Version"},
            {"name": "countSiloReads", "type": "string", "displayName": "Read Count"},
        ],
        "submissionDataTypes": {
            "files": {
                "enabled": True,
                "categories": [
                    {"name": "nucleotideAlignment"},
                    {"name": "siloReads"},
                ],
            }
        },
    }
}


def main():
    """Demonstrate schema infrastructure usage with mocked data."""
    print("=" * 60)
    print("Schema Infrastructure Demo (Mock Mode)")
    print("=" * 60)
    
    # Create schemas manually
    print("\n1. Creating SILO schema...")
    silo_schema = SiloSchema("https://lapis.example.org", "covid")
    silo_schema._schema = MOCK_SILO_SCHEMA
    
    print("   SILO metadata fields:")
    for field in silo_schema.get_metadata_fields():
        print(f"     - {field}")
    
    print("\n   SILO genes:")
    for gene in silo_schema.get_genes():
        print(f"     - {gene}")
    
    print("\n2. Creating Loculus schema...")
    loculus_schema = LoculusSchema.from_dict("covid", MOCK_LOCULUS_CONFIG)
    
    print("   Loculus metadata fields:")
    for field in loculus_schema.get_metadata_fields():
        print(f"     - {field}")
    
    print("\n3. Creating field mapping...")
    field_mapping = FieldMapping("covid")
    
    print("   Field mappings:")
    for silo_field, loculus_field in zip(
        field_mapping.get_silo_fields(),
        field_mapping.get_loculus_fields()
    ):
        print(f"     {silo_field:20s} -> {loculus_field}")
    
    print("\n4. Creating organism configuration...")
    config = OrganismConfig(
        organism="covid",
        silo_schema=silo_schema,
        loculus_schema=loculus_schema,
        field_mapping=field_mapping,
    )
    
    print("\n5. Validating configuration...")
    validation = config.validate(strict=False)
    
    if validation["valid"]:
        print("   ✓ Configuration is valid!")
    else:
        print("   ✗ Configuration has errors:")
        for error in validation["errors"]:
            print(f"     - {error}")
    
    if validation["warnings"]:
        print("\n   Warnings:")
        for warning in validation["warnings"]:
            print(f"     - {warning}")
    
    # Example: Transform SILO metadata to Loculus format
    print("\n" + "=" * 60)
    print("6. Metadata Transformation Example")
    print("=" * 60)
    
    silo_metadata = {
        "read_id": "read_12345",
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "batch001",
        "location_code": "05",
        "location_name": "Lugano (TI)",
        "sampling_date": "2024-10-08",
        "sr2silo_version": "1.5.0",
    }
    
    print("\n   SILO metadata (snake_case):")
    for key, value in silo_metadata.items():
        print(f"     {key:20s} = {value}")
    
    loculus_metadata = config.get_metadata_for_submission(silo_metadata)
    
    print("\n   Loculus metadata (camelCase):")
    for key, value in loculus_metadata.items():
        print(f"     {key:20s} = {value}")
    
    print("\n   Transformation notes:")
    print("     - Field names converted from snake_case to camelCase")
    print("     - Only mapped fields are included")
    print("     - read_id excluded (not in Loculus schema)")
    
    print("\n" + "=" * 60)
    print("Schema infrastructure demonstration complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
