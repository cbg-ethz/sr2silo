# Read Processing

Before reads end up in SILO, the Snakemake workflow subsamples each BAM file, and `sr2silo process-from-vpipe` merges paired-end reads.

## Subsampling

Subsampling runs in the `subsample_bam` rule of the workflow and is configured in the workflow config:

```yaml
ENABLE_SUBSAMPLING: true
# SUBSAMPLE_FRACTION: 0.5     # keep this fraction of records
SUBSAMPLE_MAX_READS: 4500000  # keep at most this many records
```

Depending on which keys are set, there are three modes:

1. **Fraction only:** keep `SUBSAMPLE_FRACTION` of the records.
2. **Max reads only:** keep everything if the file has at most `SUBSAMPLE_MAX_READS` records, otherwise subsample down to that number.
3. **Both:** apply the fraction, but cap the result at `SUBSAMPLE_MAX_READS`.

All current deployments (COVID, RSV-A, RSV-B) use mode 2 with a cap of 4.5 million.

Things to be aware of:

- **The cap counts BAM records, not fragments.** Records are counted with `samtools view -c`, so both mates of a pair count separately (as do secondary and supplementary alignments).
- **The cap is approximate.** Subsampling uses `samtools view -s` with a fixed seed. It keeps a random fraction of read names, so the result is close to the cap but not exact. Mates are always kept or dropped together, and reruns give the same result.
- **Logs:** per-sample details are written to `LOG_DIR/sr2silo/subsample/`, and a summary of original and subsampled counts to `LOG_DIR/subsampling/read_counts.tsv`.

## Merging paired-end reads

In paired-end sequencing, each DNA fragment is read from both ends, giving two reads (mates) with the same read name. `sr2silo process-from-vpipe` merges each pair into one read that spans the whole fragment:

1. Sort the BAM by read name, so mates are next to each other.
2. Run `paired_end_read_merger` from [smallgenomeutilities](https://github.com/cbg-ethz/smallgenomeutilities) against the nucleotide reference. Each pair becomes one record; reads without a mate are kept as they are.
3. Translate, align and write the merged reads to the output `.ndjson.zst`, one line per read.

Merging can be skipped with `--skip-merge`, but the workflow does not do this.

## Read counts

Because of merging, the read count shown in Loculus ("srSILO Read Count", the number of lines in the `.ndjson.zst`) is roughly half the number of subsampled BAM records. With a cap of 4.5 million records, the largest samples end up with about 2.25 million SILO reads.
