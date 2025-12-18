# Multi-Organism Support

sr2silo supports processing samples from multiple organisms with organism-specific reference sequences.

## Supported Organisms

| Organism | Identifier | Description |
|----------|------------|-------------|
| COVID-19 | `covid` | SARS-CoV-2 / Severe acute respiratory syndrome coronavirus 2 |
| RSV-A | `rsva` | Respiratory Syncytial Virus A |

## Reference Resolution

sr2silo resolves reference sequences in the following priority order:

1. **Local References** (fastest): `resources/references/{organism}/`
2. **LAPIS Instance**: Fetched from specified `--lapis-url` if available
3. **Fallback**: Local references as final fallback if LAPIS fetch fails

This allows for both static local references and dynamic references from a LAPIS instance.

## Usage

For detailed usage, see:
```bash
sr2silo process-from-vpipe --help
```

The `--organism` parameter specifies which organism to process. Can also be set via `ORGANISM` environment variable (CLI argument takes precedence).

## Adding New Organisms

To add support for a new organism:

1. **Add reference files** to `resources/references/{organism_id}/`:
   - `nuc_ref.fasta` - Nucleotide reference sequence(s)
   - `aa_ref.fasta` - Amino acid reference sequences for gene annotations

2. **Use GenBank parser** (if starting from GenBank format):
   ```bash
   python scripts/extract_gbk_references.py \
       reference.gbk \
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

Specify organism in `workflow/config.yaml`:
```yaml
ORGANISM: "covid"  # or "rsva"
```

## Troubleshooting

**Reference files not found:**
- Verify files exist: `ls resources/references/{organism}/`
- Check organism identifier spelling
- Use `--lapis-url` to fetch from LAPIS
- Generate from GenBank: `python scripts/extract_gbk_references.py --help`
