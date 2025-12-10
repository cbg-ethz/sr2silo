# Multi-Organism Support

sr2silo supports processing samples from multiple organisms with organism-specific reference sequences.

## Supported Organisms

| Organism | Identifier | Description |
|----------|------------|-------------|
| SARS-CoV-2 | `sars-cov-2` | Severe acute respiratory syndrome coronavirus 2 (COVID-19) |
| RSV-A | `rsva` | Respiratory Syncytial Virus A |

## Reference Resolution

sr2silo resolves reference sequences in the following priority order:

1. **Local References** (fastest): `resources/references/{organism}/`
2. **LAPIS Instance**: Fetched from specified `--lapis-url` if available
3. **Fallback**: Local references as final fallback if LAPIS fetch fails

This allows for both static local references and dynamic references from a LAPIS instance.

## Usage

### Specifying Organism via CLI

```bash
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --organism sars-cov-2 \
    --output-fp output.ndjson.zst
```

### Using Environment Variable

Set the `ORGANISM` environment variable instead of passing `--organism`:

```bash
export ORGANISM=rsva

sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --output-fp output.ndjson.zst
```

**Note:** CLI `--organism` parameter overrides the `ORGANISM` environment variable.

### With LAPIS Reference Fetching

To fetch references from a LAPIS instance (with local fallback):

```bash
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --organism rsva \
    --lapis-url https://lapis.example.com \
    --output-fp output.ndjson.zst
```

The system will:
1. Check for local RSV-A references first
2. If not found, try fetching from the LAPIS instance
3. Use local references as fallback if LAPIS fetch fails

## Adding New Organisms

To add support for a new organism:

1. **Add reference files** to `resources/references/{organism_id}/`:
   - `nuc_ref.fasta` - Nucleotide reference sequence(s)
   - `aa_ref.fasta` - Amino acid reference sequences for gene annotations

2. **Use GenBank parser** (if starting from GenBank format):
   ```bash
   python scripts/extract_gbk_references.py \
       --gbk-file reference.gbk \
       --nuc-output resources/references/{organism_id}/nuc_ref.fasta \
       --aa-output resources/references/{organism_id}/aa_ref.fasta
   ```

3. **Update configuration** (optional):
   - Update workflow `config.yaml` with your new organism ID
   - Update documentation with organism details

4. **Test** (optional):
   - Add test fixtures in `tests/conftest.py`
   - Add parameterized tests in `tests/test_main.py`

## Workflow Integration

When using the Snakemake workflow, specify organism in `workflow/config.yaml`:

```yaml
# Organism identifier for processing
ORGANISM: "sars-cov-2"  # or "rsva", etc.
```

Or override at runtime:

```bash
snakemake --config ORGANISM=rsva
```

## Troubleshooting

### Reference Files Not Found

If you get an error about missing reference files:

```
FileNotFoundError: No reference files found for organism 'organism_id'
```

**Solutions:**
1. Verify reference files exist: `ls resources/references/{organism}/`
2. Check organism name spelling matches exactly
3. Use `--lapis-url` to fetch references from LAPIS if local files are unavailable
4. Use scripts/extract_gbk_references.py to generate references from GenBank files

### Mismatched References

If your processing results seem incorrect:

1. Verify organism identifier matches your data source
2. Check reference file versions are appropriate for your samples
3. Consider using `--lapis-url` to ensure you have the latest references
