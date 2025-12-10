# Configuration Guide

sr2silo supports flexible configuration through command-line arguments and environment variables, making it easy to use in different deployment scenarios.

## Command-Line Arguments

View all available arguments:

```bash
sr2silo process-from-vpipe --help
sr2silo submit-to-loculus --help
```

## Environment Variables

Environment variables provide configuration defaults. CLI arguments override environment variables.

### Processing Configuration (process-from-vpipe)

| Variable | Purpose | Default | Example |
|----------|---------|---------|---------|
| `ORGANISM` | Organism identifier | None (required) | `sars-cov-2`, `rsva` |
| `TIMELINE_FILE` | Metadata timeline file | None (required) | `/path/to/timeline.tsv` |
| `LAPIS_URL` | LAPIS instance URL (optional) | None | `https://lapis.example.com` |

### Submission Configuration (submit-to-loculus)

| Variable | Purpose | Default | Example |
|----------|---------|---------|---------|
| `ORGANISM` | Organism identifier | None (required) | `sars-cov-2`, `rsva` |
| `KEYCLOAK_TOKEN_URL` | Authentication endpoint | None (required) | `https://auth.example.com/token` |
| `BACKEND_URL` | SILO backend API endpoint | None (required) | `https://api.example.com/api` |
| `GROUP_ID` | Loculus group ID | None (required) | `1`, `42` |
| `USERNAME` | Submission username | None (required) | Your username |
| `PASSWORD` | Submission password | None (required) | Your password |

## Usage Examples

### Example 1: Local Processing with Environment Variables

```bash
# Set environment variables
export ORGANISM=sars-cov-2
export TIMELINE_FILE=/path/to/timeline.tsv

# Process data (uses environment variables)
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --output-fp output.ndjson.zst
```

### Example 2: CLI Arguments Override Environment

```bash
# Set default organism
export ORGANISM=sars-cov-2

# Override with CLI argument for RSV data
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file /path/to/rsv_timeline.tsv \
    --organism rsva \
    --output-fp output.ndjson.zst
```

### Example 3: Complete Processing and Submission Workflow

```bash
# Set up environment
export ORGANISM=sars-cov-2
export KEYCLOAK_TOKEN_URL=https://auth.db.wasap.genspectrum.org/realms/loculus/protocol/openid-connect/token
export BACKEND_URL=https://api.db.wasap.genspectrum.org/backend
export GROUP_ID=1
export USERNAME=your-username
export PASSWORD=your-password

# Step 1: Process data
sr2silo process-from-vpipe \
    --input-file aligned_reads.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --output-fp sample_001.ndjson.zst

# Step 2: Submit to Loculus
sr2silo submit-to-loculus \
    --processed-file sample_001.ndjson.zst \
    --nucleotide-alignment aligned_reads.bam
```

### Example 4: Using LAPIS for Dynamic References

```bash
# Process with LAPIS for reference fetching
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --organism rsva \
    --lapis-url https://lapis.example.com \
    --output-fp output.ndjson.zst
```

## Configuration in Workflows

### Snakemake Workflow

Configure in `workflow/config.yaml`:

```yaml
# Organism identifier
ORGANISM: "sars-cov-2"

# LAPIS configuration (optional)
LAPIS_URL: "https://lapis.wasap.genspectrum.org/"

# Submission credentials
KEYCLOAK_TOKEN_URL: "https://auth.db.wasap.genspectrum.org/..."
BACKEND_URL: "https://api.db.wasap.genspectrum.org/..."
GROUP_ID: 1
```

Or override at runtime:

```bash
snakemake --config ORGANISM=rsva LAPIS_URL="https://custom.lapis.com"
```

## Best Practices

1. **Use environment variables for secrets**: Store credentials in environment variables, not in code or configuration files
   ```bash
   export KEYCLOAK_TOKEN_URL="..."
   export USERNAME="..."
   export PASSWORD="..."
   ```

2. **Use CLI arguments for experiment-specific settings**: Override defaults for specific runs
   ```bash
   sr2silo process-from-vpipe \
       --organism rsva \  # Override default
       --lapis-url https://custom.lapis.com  # Use custom LAPIS
   ```

3. **Document your environment setup**: Create a setup script for your deployment
   ```bash
   # setup_env.sh
   export ORGANISM=sars-cov-2
   export KEYCLOAK_TOKEN_URL=https://auth.example.com/token
   export BACKEND_URL=https://api.example.com/api
   source setup_env.sh
   ```

4. **Verify configuration before processing**: Use `--help` to check what values will be used
   ```bash
   sr2silo process-from-vpipe --help
   ```

## Troubleshooting

### "ORGANISM environment variable is not set"

**Error:** Required `ORGANISM` variable missing

**Solutions:**
- Set via environment: `export ORGANISM=sars-cov-2`
- Set via CLI: `--organism sars-cov-2`

### "Reference files not found"

**Error:** sr2silo can't locate reference sequences

**Solutions:**
- Check organism identifier spelling
- Verify reference files exist: `ls resources/references/{organism}/`
- Use `--lapis-url` to fetch from LAPIS
- See [Multi-Organism Support](organisms.md) for adding new organisms

### "Authentication failed"

**Error:** Submission to Loculus fails with credentials

**Solutions:**
- Verify credentials in environment variables
- Check endpoint URLs are correct
- Ensure GROUP_ID is valid for your account
