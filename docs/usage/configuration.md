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

## Snakemake Workflow

Configure in `workflow/config.yaml`:
```yaml
ORGANISM: "sars-cov-2"
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
