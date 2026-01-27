# Data Schemas

Pydantic schemas defining the SILO database format for reads.

These schemas validate read data before submission to SILO, ensuring all records conform to the expected format. Sequences follow [IUPAC nucleotide codes](https://www.bioinformatics.org/sms2/iupac.html).

## Output Format

The output format for LAPIS-SILO v0.8.0+ uses:

- **camelCase naming** for metadata fields (e.g., `readId`, `sampleId`, `batchId`)
- **Flat structure** with metadata at root level (no nested "metadata" object)
- **Segment format** with `sequence`, `insertions`, and `offset` fields
- **Insertion format** as `"position:sequence"` (e.g., `"123:ACGT"`)

## Modifying the Schema

To customize the output schema for your use case:

1. **Edit the schema file** at `src/sr2silo/silo_read_schema.py`:
   - Add or modify fields in the `ReadMetadata` class
   - Use Pydantic `Field(alias="camelCaseName")` for camelCase output

2. **Update the database config** at `resources/silo/database_config.yaml`:
   - Ensure field names match the Pydantic aliases

3. **Run validation** to verify consistency:
   ```bash
   python tests/test_database_config_validation.py
   ```

This validation ensures your Pydantic schema matches the SILO database configuration.

## ReadMetadata

::: sr2silo.silo_read_schema.ReadMetadata

## NucleotideSegment

::: sr2silo.silo_read_schema.NucleotideSegment

## AminoAcidSegment

::: sr2silo.silo_read_schema.AminoAcidSegment

## AlignedReadSchema

::: sr2silo.silo_read_schema.AlignedReadSchema
