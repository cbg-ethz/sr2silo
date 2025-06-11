# Environment Variable Configuration

sr2silo supports configuration via environment variables, making it easy for users who install via conda or pip to set up default values without needing to specify all parameters on every command line invocation.

## How It Works

Environment variables serve as **defaults** for command-line arguments. The precedence is:

1. **CLI arguments** (highest priority) - always override environment variables
2. **Environment variables** - provide defaults when CLI arguments are not specified  
3. **Built-in defaults** (lowest priority) - used when neither CLI args nor env vars are set

## Environment Variables

### For `process-from-vpipe` command:

| Environment Variable | CLI Argument | Description | Example |
|---------------------|--------------|-------------|---------|
| `SAMPLE_DIR` | `--input-file` | Directory containing `REF_aln_trim.bam` | `./data/sample/alignments/` |
| `SAMPLE_ID` | `--sample-id` | Sample identifier | `A1_05_2024_10_08` |
| `BATCH_ID` | `--batch-id` | Batch identifier | `20241024_2411515907` |
| `TIMELINE_FILE` | `--timeline-file` | Path to timeline metadata file | `./data/timeline.tsv` |
| `PRIMER_FILE` | `--primer-file` | Path to primers configuration file | `./data/primers.yaml` |
| `RESULTS_DIR` | `--output-fp` | Results directory (auto-generates filename) | `./results/` |
| `NEXTCLADE_REFERENCE` | `--reference` | Reference genome (default: `sars-cov-2`) | `sars-cov-2` |

### For `submit-to-loculus` command:

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `KEYCLOAK_TOKEN_URL` | Authentication endpoint | `https://authentication-wise-seqs.loculus.org/...` |
| `SUBMISSION_URL` | SILO submission endpoint | `https://backend-wise-seqs.loculus.org/...` |

### System Variables:

| Environment Variable | Description | Default |
|---------------------|-------------|---------|
| `CI` | Set to `true` in CI environments | `false` |

## Usage Examples

### Setting up environment variables

You can set environment variables in several ways:

**1. In your shell profile** (persistent):
```bash
# Add to ~/.bashrc, ~/.zshrc, etc.
export TIMELINE_FILE="$HOME/sr2silo/timeline.tsv"
export PRIMER_FILE="$HOME/sr2silo/primers.yaml"
export NEXTCLADE_REFERENCE="sars-cov-2"
export KEYCLOAK_TOKEN_URL="https://your-auth-server.com/token"
export SUBMISSION_URL="https://your-silo-server.com/submit"
```

**2. Using a .env file** (in your working directory):
```bash
# Create .env file in your working directory
cat > .env << EOF
TIMELINE_FILE=./data/timeline.tsv
PRIMER_FILE=./data/primers.yaml
NEXTCLADE_REFERENCE=sars-cov-2
KEYCLOAK_TOKEN_URL=https://your-auth-server.com/token
SUBMISSION_URL=https://your-silo-server.com/submit
EOF

# Load environment variables (bash/zsh)
set -a; source .env; set +a
```

**3. Per-command** (temporary):
```bash
SAMPLE_DIR="/data/sample1" SAMPLE_ID="A1_05_2024_10_08" sr2silo process-from-vpipe
```

### Running commands with environment variables

**With environment variables providing defaults:**
```bash
# Set common configuration
export TIMELINE_FILE="./data/timeline.tsv"
export PRIMER_FILE="./data/primers.yaml" 
export RESULTS_DIR="./results/"

# Now you only need to specify per-run parameters
sr2silo process-from-vpipe \
    --input-file ./data/sample1/alignments/REF_aln_trim.bam \
    --sample-id A1_05_2024_10_08 \
    --batch-id 20241024_2411515907
```

**Using SAMPLE_DIR for automatic input file detection:**
```bash
# Set the sample directory - sr2silo will look for REF_aln_trim.bam inside it
export SAMPLE_DIR="./data/sample1/alignments/"
export SAMPLE_ID="A1_05_2024_10_08"
export BATCH_ID="20241024_2411515907"
export TIMELINE_FILE="./data/timeline.tsv"
export PRIMER_FILE="./data/primers.yaml"
export RESULTS_DIR="./results/"

# Now you can run with no arguments!
sr2silo process-from-vpipe
```

**CLI arguments override environment variables:**
```bash
# Even if SAMPLE_ID is set in environment, CLI argument takes precedence
export SAMPLE_ID="env_sample"
sr2silo process-from-vpipe --sample-id "cli_sample" # Uses "cli_sample"
```

## Installation-specific Setup

### For conda users:
After installing via conda, create a configuration script:
```bash
# Create ~/sr2silo-config.sh
cat > ~/sr2silo-config.sh << 'EOF'
#!/bin/bash
export TIMELINE_FILE="$HOME/sr2silo-data/timeline.tsv"
export PRIMER_FILE="$HOME/sr2silo-data/primers.yaml"
export NEXTCLADE_REFERENCE="sars-cov-2"
export KEYCLOAK_TOKEN_URL="https://your-server.com/token"
export SUBMISSION_URL="https://your-server.com/submit"
EOF

# Source it when needed
source ~/sr2silo-config.sh
sr2silo process-from-vpipe --input-file ./sample.bam # other args from env
```

### For pip users:
Similar to conda, but you might want to integrate with your Python virtual environment:
```bash
# Add to your virtual environment's activate script
echo 'export TIMELINE_FILE="..."' >> $VIRTUAL_ENV/bin/activate
```

## Best Practices

1. **Set installation-time variables once**: `TIMELINE_FILE`, `PRIMER_FILE`, `NEXTCLADE_REFERENCE`, `KEYCLOAK_TOKEN_URL`, `SUBMISSION_URL`

2. **Set per-run variables as needed**: `SAMPLE_DIR`, `SAMPLE_ID`, `BATCH_ID`, `RESULTS_DIR`

3. **Use CLI arguments for one-off overrides**: When you need to override a default for a single run

4. **Use .env files for project-specific settings**: Keep different configurations for different projects

5. **Check your configuration**: Use `sr2silo process-from-vpipe --help` to see which environment variables are supported