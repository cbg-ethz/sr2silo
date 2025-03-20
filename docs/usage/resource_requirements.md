# Resource Requirements

When running sr2silo, especially the `--import to loculus` command, you should be aware of the following resource requirements:

## Memory Requirements

- **Standard Resources**: The default Snakemake configuration uses 8 GB RAM and one CPU core
- **Actual Usage**: 
  - sr2silo processes in batches of 100k reads at a time, requiring approximately 3 GB of RAM
  - An additional 3 GB of RAM is needed for running Diamond for amino acid translation and alignment

## Storage Requirements

- **Temporary Space**: 
  - sr2silo itself requires some temporary directory space
  - Diamond may require up to 30 GB of temporary storage
  - Ideally, this should be on high I/O storage

## Cluster Environment Configuration

When running on a personal computer, standard temporary directories are usually sufficient as long as you have free disk space. However, in a cluster environment:

1. Set the environment variable `TMPDIR` to a location with at least 50 GB of free space per run:

```bash
export TMPDIR=/path/to/temp/directory
```

2. Make sure this directory has good I/O performance for optimal processing speed

This configuration will help prevent issues with temporary file storage during the processing of large datasets.
