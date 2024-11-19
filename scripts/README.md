# Docker: Process a Single BAM File

## Configuration

Edit the `docker-compose.env` file in the `docker-compose` directory with the following paths:

```env
SAMPLE_DIR=../../../data/sr2silo/daemon_test/samples/A1_05_2024_10_08/20241024_2411515907/alignments/
SAMPLE_ID=A1_05_2024_10_08
BATCH_ID=20241024_2411515907
TIMELINE_FILE=../../../data/sr2silo/daemon_test/timeline.tsv
NEXTCLADE_REFERENCE=sars-cov2
RESULTS_DIR=./results
```

## Processing

To process a single sample, run the following command:

```sh
docker-compose --env-file docker-compose.env up --build
```
