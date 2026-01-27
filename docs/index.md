# sr2silo

**Convert BAM nucleotide alignments to cleartext alignments for LAPIS-SILO**

sr2silo processes short-read nucleotide alignments from `.bam` files, translates and aligns reads in amino acids, and outputs JSON compatible with [LAPIS-SILO](https://github.com/GenSpectrum/LAPIS-SILO) v0.8.0+.

## Features

- **BAM to JSON conversion**: Transform CIGAR alignments to cleartext format
- **Amino acid translation**: Translate and align reads using Diamond/BlastX
- **Insertion/deletion handling**: Gracefully extract and represent indels
- **Multi-organism support**: Process COVID-19, RSV-A, and other pathogens
- **Loculus integration**: Direct submission to Loculus backend

## Quick Start

```bash
# Install from Bioconda
conda install -c bioconda sr2silo

# Process BAM data
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --organism covid \
    --output-fp output.ndjson.zst

# Submit to Loculus
sr2silo submit-to-loculus \
    --processed-file output.ndjson.zst
```

## Input/Output

**Input**: BAM/SAM file with CIGAR alignments:
```text
294 163 NC_045512.2 79  60  31S220M =   197 400 CTCTTGTAGAT FGGGHHHHLMM ...
```

**Output**: Newline-delimited JSON (`.ndjson.zst`) compatible with LAPIS-SILO v0.8.0+:

```json
{
  "readId": "AV233803:AV044:2411515907:1:10805:5199:3294",
  "sampleId": "A1_05_2024_10_08",
  "batchId": "20241024_2411515907",
  "samplingDate": "2024-10-08",
  "locationName": "Lugano (TI)",
  "locationCode": "5",
  "sr2siloVersion": "1.3.0",
  "main": {
    "sequence": "CGGTTTCGTCCGTGTTGCAGCCG...GTGTCAACATCTTAAAGATGGCACTTGTG",
    "insertions": ["10:ACTG", "456:TACG"],
    "offset": 4545
  },
  "S": {
    "sequence": "MESLVPGFNEKTHVQLSLPVLQVRVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGV",
    "insertions": ["23:A", "145:KLM"],
    "offset": 78
  },
  "ORF1a": {
    "sequence": "XXXMESLVPGFNEKTHVQLSLPVLQVRVRGFGDSVEEVLSEARQHLKDGTCGLV",
    "insertions": ["2323:TG", "2389:CA"],
    "offset": 678
  },
  "E": null,
  "M": null,
  "N": null
}
```

See [Data Schemas](api/schema.md) for details on the output format and how to customize it.

## Documentation

- [Configuration](usage/configuration.md) - Environment variables and CLI options
- [Multi-Organism Support](usage/organisms.md) - Supported organisms and adding new ones
- [Deployment](usage/deployment.md) - Multi-virus cluster deployment
- [Resource Requirements](usage/resource_requirements.md) - Memory and storage needs

## Links

- [GitHub Repository](https://github.com/cbg-ethz/sr2silo)
- [Bioconda Package](https://bioconda.github.io/recipes/sr2silo/README.html)
