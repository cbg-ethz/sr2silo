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
| `ORGANISM` | Organism identifier | None (required) | `covid`, `rsva` |
| `TIMELINE_FILE` | Metadata timeline file | None (required) | `/path/to/timeline.tsv` |
| `LAPIS_URL` | LAPIS instance URL (optional) | None | `https://lapis.example.com` |
| `REFERENCE_ACCESSION` | Filter reads by reference accession | `""` (all reads) | `EPI_ISL_412866` |

#### Reference Filtering

When processing BAM files containing alignments to multiple similar references (e.g., RSV-A and RSV-B), use `--reference-accession` to filter reads:

```bash
sr2silo process-from-vpipe \
  --input-file alignments.bam \
  --reference-accession "EPI_ISL_412866" \
  ...
```

To find available reference accessions in your BAM file:

```bash
samtools view -H file.bam | grep @SQ
```

If filtering results in zero reads, processing terminates with a `ZeroFilteredReadsError`.

### Submission Configuration (submit-to-loculus)

| Variable | Purpose | Default | Example |
|----------|---------|---------|---------|
| `ORGANISM` | Organism identifier | None (required) | `covid`, `rsva` |
| `KEYCLOAK_TOKEN_URL` | Authentication endpoint | None (required) | `https://auth.example.com/token` |
| `BACKEND_URL` | SILO backend API endpoint | None (required) | `https://api.example.com/api` |
| `GROUP_ID` | Loculus group ID | None (required) | `1`, `42` |
| `USERNAME` | Submission username | None (required) | Your username |
| `PASSWORD` | Submission password | None (required) | Your password |
| `AUTO_RELEASE` | Auto-approve sequences after submission | `false` | `true` |

#### Auto-Release

Automatically approve/release sequences after submission:

```bash
# Via CLI flag
sr2silo submit-to-loculus \
  -f data.ndjson.zst \
  -a alignment.bam \
  --auto-release

# With custom delay (default: 180 seconds)
sr2silo submit-to-loculus \
  -f data.ndjson.zst \
  -a alignment.bam \
  --auto-release \
  --release-delay 60

# Via environment variable
export AUTO_RELEASE=true
sr2silo submit-to-loculus -f data.ndjson.zst -a alignment.bam
```

The `--release-delay` option (default 180s) allows backend processing time before release.

## Snakemake Workflow

Configure in `workflow/config.yaml`:
```yaml
ORGANISM: "covid"
REFERENCE_ACCESSION: ""  # Filter reads by reference (use for RSV-A/B)
AUTO_RELEASE: false      # Auto-approve after submission
LAPIS_URL: "https://lapis.wasap.genspectrum.org/"
KEYCLOAK_TOKEN_URL: "https://auth.db.wasap.genspectrum.org/..."
BACKEND_URL: "https://api.db.wasap.genspectrum.org/..."
GROUP_ID: 1
```

## Best Practices

- **Store credentials in environment variables**, not in code
- **Use CLI arguments for experiment-specific settings** that override defaults
- **Create setup scripts** for your deployment environment
- **Use `--help`** to verify available options

## Troubleshooting

**ORGANISM not set:** Export variable or use `--organism` flag

**Reference files not found:** Check spelling, verify files exist, or use `--lapis-url`

**Authentication failed:** Verify credentials and endpoint URLs in environment variables
