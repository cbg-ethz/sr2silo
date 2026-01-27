# Reference File Examples

This directory contains example reference files for sr2silo.
These files are **NOT** included in the installed package.

## Usage

### Option 1: Explicit Paths (Recommended for Production)

Provide explicit paths to reference files:

```bash
sr2silo process-from-vpipe \
  --nuc-ref examples/references/covid/nuc_ref.fasta \
  --aa-ref examples/references/covid/aa_ref.fasta \
  --input-file alignments.bam \
  --sample-id sample123 \
  --timeline-file timeline.tsv \
  --output-fp output.ndjson.zst
```

### Option 2: LAPIS Auto-Fetch

Fetch references automatically from a LAPIS instance:

```bash
sr2silo process-from-vpipe \
  --lapis-url https://lapis.wasap.genspectrum.org/covid \
  --input-file alignments.bam \
  --sample-id sample123 \
  --timeline-file timeline.tsv \
  --output-fp output.ndjson.zst
```

References are cached to the Python package's resources directory
(e.g., `<site-packages>/resources/references/<lapis-hostname>/`) rather than `~/.cache/sr2silo/references/`.
When using conda environments, this will be inside the environment's `lib/pythonX.XX/resources/references/` folder.

### Option 3: Environment Variables

```bash
export NUC_REF=/path/to/nuc_ref.fasta
export AA_REF=/path/to/aa_ref.fasta
sr2silo process-from-vpipe ...
```

## File Format

Each organism directory should contain:

- `nuc_ref.fasta`: Nucleotide reference sequence (FASTA format)
- `aa_ref.fasta`: Amino acid reference sequences for genes (FASTA format)

### Example Structure

```
references/
├── covid/
│   ├── nuc_ref.fasta    # SARS-CoV-2 genome (NC_045512.2)
│   └── aa_ref.fasta     # Protein sequences for genes
└── rsva/
    ├── nuc_ref.fasta    # RSV-A genome
    └── aa_ref.fasta     # Protein sequences for genes
```

## Priority Order

When resolving references, sr2silo uses this priority:

1. CLI flags (`--nuc-ref`, `--aa-ref`)
2. Environment variables (`NUC_REF`, `AA_REF`)
3. LAPIS auto-fetch with cache (if `--lapis-url` provided)
