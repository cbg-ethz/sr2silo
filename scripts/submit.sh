#!/bin/bash

sbatch \
    --mail-type=END \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=1600 \
    --time=2:00:00 \
    -o /cluster/home/koehng/logs/sr2silo.out \
    -e /cluster/home/koehng/logs/sr2silo.err \
    ./run.sh
