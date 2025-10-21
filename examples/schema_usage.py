"""Example usage of the schema infrastructure.

This example demonstrates how to:
1. Fetch SILO schema from API with caching
2. Load Loculus schema from configuration
3. Set up organism-specific field mapping
4. Validate the complete configuration
5. Transform metadata for submission
"""

from sr2silo.schema import OrganismConfig

# Example Loculus configuration (similar to WisePulse)
LOCULUS_CONFIG = {
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
    """Demonstrate schema infrastructure usage."""
    # Create organism configuration by fetching SILO schema and loading Loculus config
    print("Setting up organism configuration for covid...")
    
    config = OrganismConfig.from_api_and_config(
        organism="covid",
        lapis_url="https://lapis.wasap.genspectrum.org",
        loculus_config=LOCULUS_CONFIG,
        use_cache=True,  # Use cached schema if available
    )
    
    # Log configuration summary
    print("\n" + "=" * 60)
    config.log_configuration_summary()
    print("=" * 60 + "\n")
    
    # Validate the configuration
    print("Validating configuration...")
    validation = config.validate(strict=False)
    
    if validation["valid"]:
        print("✓ Configuration is valid!")
    else:
        print("✗ Configuration has errors:")
        for error in validation["errors"]:
            print(f"  - {error}")
    
    if validation["warnings"]:
        print("\nWarnings:")
        for warning in validation["warnings"]:
            print(f"  - {warning}")
    
    # Example: Transform SILO metadata to Loculus format
    print("\n" + "=" * 60)
    print("Example metadata transformation:")
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
    
    print("\nSILO metadata (snake_case):")
    for key, value in silo_metadata.items():
        print(f"  {key}: {value}")
    
    loculus_metadata = config.get_metadata_for_submission(silo_metadata)
    
    print("\nLoculus metadata (camelCase):")
    for key, value in loculus_metadata.items():
        print(f"  {key}: {value}")
    
    print("\n" + "=" * 60)
    print("Schema infrastructure demonstration complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
