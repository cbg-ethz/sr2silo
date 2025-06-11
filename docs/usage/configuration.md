# Configuration

sr2silo supports configuration through both command-line arguments and environment variables. Environment variables provide a convenient way to set defaults for parameters that don't change frequently between runs.

## Environment Variables

The following environment variables can be used to configure sr2silo:

### Processing Configuration

- `TIMELINE_FILE`: Path to the timeline file (used by `process-from-vpipe`)
- `PRIMER_FILE`: Path to the primers file (used by `process-from-vpipe`)  
- `NEXTCLADE_REFERENCE`: Reference genome identifier (used by `process-from-vpipe`)

### Submission Configuration

- `KEYCLOAK_TOKEN_URL`: Keycloak authentication URL (used by `submit-to-loculus`)
- `SUBMISSION_URL`: Loculus submission URL (used by `submit-to-loculus`)

### System Configuration

- `CI`: Set to "true" to enable CI environment mode
- `TMPDIR`: Temporary directory for intermediate files

## Usage Patterns

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

You can also use a `.env` file to set environment variables:

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
sr2silo process-from-vpipe --input-file /path/to/input.bam --sample-id SAMPLE_001 --batch-id BATCH_001 --output-fp /path/to/output.ndjson
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

## Parameter Categories

### Run-specific Parameters
These parameters typically change for each run and should be provided via command-line arguments:
- `--input-file` / `-i`: Input BAM file
- `--sample-id` / `-s`: Sample identifier  
- `--batch-id` / `-b`: Batch identifier
- `--output-fp` / `-o`: Output file path

### Installation-specific Parameters  
These parameters are typically set once after installation and can be configured via environment variables:
- `--timeline-file` / `-t`: Timeline file path
- `--primer-file` / `-p`: Primers file path
- `--reference` / `-r`: Reference genome

### System Parameters
These parameters control system behavior and are usually set via environment variables:
- `TMPDIR`: Temporary directory
- `CI`: CI environment flag