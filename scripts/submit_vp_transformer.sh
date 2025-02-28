#!/bin/bash

sbatch \
    --mail-type=END \
    --ntasks=1 \
    --cpus-per-task=8 \
    --mem-per-cpu=8000 \
    --time=6:00:00 \
    -o /cluster/home/koehng/logs/sr2silo.out \
    -e /cluster/home/koehng/logs/sr2silo.err \
    ./run_vp_transformer.sh
