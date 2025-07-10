# Environment Variable Configuration

sr2silo supports flexible configuration through environment variables with CLI parameter override capability. This is especially useful for users who install sr2silo via conda packages or as a Python package where managing a `.env` file is not practical.

## Overview

The package works without requiring a `.env` file. Configuration can be provided through:

1. **Environment variables** - Set in your shell or system environment
2. **CLI parameters** - Provided directly to commands
3. **Built-in defaults** - Fallback values where applicable

**Precedence order**: CLI parameters > Environment variables > Built-in defaults

## Configuration Parameters

### For `process-from-vpipe` command

| CLI Parameter | Environment Variable | Description | Default |
|---------------|---------------------|-------------|---------|
| `--timeline-file` | `TIMELINE_FILE` | Path to timeline file | *(required)* |
| `--primer-file` | `PRIMER_FILE` | **DEPRECATED** - All metadata now sourced from timeline | *(ignored)* |
| `--reference` | `NEXTCLADE_REFERENCE` | Reference genome identifier | `sars-cov-2` |

### For `submit-to-loculus` command

| CLI Parameter | Environment Variable | Description | Default |
|---------------|---------------------|-------------|---------|
| `--keycloak-token-url` | `KEYCLOAK_TOKEN_URL` | Keycloak authentication URL | *(required)* |
| `--submission-url` | `SUBMISSION_URL` | Loculus submission URL | *(required)* |

## Usage Patterns

### Using Environment Variables

Set environment variables to provide defaults:

```bash
export TIMELINE_FILE=/path/to/timeline.tsv
export NEXTCLADE_REFERENCE=sars-cov-2

# Now you can run without specifying these parameters
sr2silo process-from-vpipe \
  --input-file /path/to/input.bam \
  --sample-id SAMPLE_001 \
  --batch-id BATCH_001 \
  --output-fp /path/to/output.ndjson
```

### Using Configuration Files

You can also use a `.env` file to set environment variables:

```bash
# .env file
TIMELINE_FILE=/path/to/timeline.tsv
NEXTCLADE_REFERENCE=sars-cov-2
KEYCLOAK_TOKEN_URL=https://your-keycloak-url/token
SUBMISSION_URL=https://your-submission-url
```

Then source the file before running sr2silo:

```bash
source .env
sr2silo process-from-vpipe \
  --input-file /path/to/input.bam \
  --sample-id SAMPLE_001 \
  --batch-id BATCH_001 \
  --output-fp /path/to/output.ndjson
```

### Overriding Environment Variables

Command-line arguments always take precedence over environment variables:

```bash
export NEXTCLADE_REFERENCE=sars-cov-2

# This will use "custom-ref" instead of "sars-cov-2"
sr2silo process-from-vpipe \
  --reference custom-ref \
  --input-file /path/to/input.bam \
  --sample-id SAMPLE_001 \
  --batch-id BATCH_001 \
  --output-fp /path/to/output.ndjson
```

### Submit to Loculus Configuration

For the submit command, configure the authentication and submission URLs:

```bash
export KEYCLOAK_TOKEN_URL=https://authentication.example.com/token
export SUBMISSION_URL=https://backend.example.com/submit

sr2silo submit-to-loculus \
  --processed-file data.ndjson.zst \
  --sample-id SAMPLE_001
```

Or override with CLI parameters:

```bash
sr2silo submit-to-loculus \
  --processed-file data.ndjson.zst \
  --sample-id SAMPLE_001 \
  --keycloak-token-url https://auth.example.com/token \
  --submission-url https://backend.example.com/submit
```