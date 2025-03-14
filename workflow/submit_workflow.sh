#!/bin/bash
sbatch \
    --mail-type=END \
    --ntasks=1 \
    --cpus-per-task=12 \
    --mem-per-cpu=16000 \
    --time=2:00:00 \
    -o /cluster/home/koehng/logs/sr2silo.out \
    -e /cluster/home/koehng/logs/sr2silo.err \
    snakemake -c 12
