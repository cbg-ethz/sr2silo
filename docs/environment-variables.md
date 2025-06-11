# Environment Variable Configuration

This document describes how to configure sr2silo using environment variables and command-line arguments.

## Overview

sr2silo now supports flexible configuration through environment variables that can be overridden by command-line arguments. The precedence order is:

1. Command-line arguments (highest priority)
2. Environment variables  
3. Built-in defaults (lowest priority)

## Environment Variables

### For `process-from-vpipe` command:

- `TIMELINE_FILE`: Path to the timeline file (TSV format)
- `PRIMER_FILE`: Path to the primers file (YAML format)  
- `NEXTCLADE_REFERENCE`: Reference genome name (default: "sars-cov-2")

### For `submit-to-loculus` command:

- `KEYCLOAK_TOKEN_URL`: Keycloak authentication URL
- `SUBMISSION_URL`: Loculus submission URL

## Usage Examples

### Using Environment Variables

Set environment variables to provide defaults:

```bash
export TIMELINE_FILE=/path/to/timeline.tsv
export PRIMER_FILE=/path/to/primers.yaml
export NEXTCLADE_REFERENCE=sars-cov-2

# Now you can run without specifying these parameters
sr2silo process-from-vpipe \
  --input-file /path/to/input.bam \
  --sample-id SAMPLE_001 \
  --batch-id BATCH_001 \
  --output-fp /path/to/output.ndjson
```

### Using Configuration Files

You can also use a .env file to set environment variables:

```bash
# .env file
TIMELINE_FILE=/path/to/timeline.tsv
PRIMER_FILE=/path/to/primers.yaml
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

### Submit to Loculus with URLs

```bash
export KEYCLOAK_TOKEN_URL=https://auth.example.com/token
export SUBMISSION_URL=https://backend.example.com/submit

sr2silo submit-to-loculus \
  --processed-file /path/to/processed.ndjson.zst \
  --sample-id SAMPLE_001
  
# Or override with CLI arguments:
sr2silo submit-to-loculus \
  --processed-file /path/to/processed.ndjson.zst \
  --sample-id SAMPLE_001 \
  --keycloak-token-url https://different-auth.example.com/token \
  --submission-url https://different-backend.example.com/submit
```

## Installation Types

### Conda Package Installation

When installing via conda, you can set environment variables in your shell profile:

```bash
echo 'export TIMELINE_FILE=/path/to/timeline.tsv' >> ~/.bashrc
echo 'export PRIMER_FILE=/path/to/primers.yaml' >> ~/.bashrc
source ~/.bashrc
```

### Python Package Installation

When installing as a Python package, you can:

1. Set environment variables globally
2. Use a .env file in your working directory
3. Pass all parameters via command-line arguments

## Migration from .env File

If you were previously using a .env file, you can:

1. Keep using it (sr2silo will still load it if present)
2. Convert to system environment variables
3. Use command-line arguments for all parameters

The package now works without requiring a .env file, making it more portable and suitable for different deployment scenarios.