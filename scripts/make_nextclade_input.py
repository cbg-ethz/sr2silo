"""Make cleartext alignment from BAM file with reference sequence
    for Nextclade input."""

from __future__ import annotations

import json
import time
from pathlib import Path
from statistics import mean, variance

import nextclade
from sr2silo.process.convert import bam_to_cleartext_alignment

REFERENCE_PATH = Path("resources/sars-cov-2/NC_045512.2.fasta")
OUTPUT_PATH = Path("out/nextclade/nextclade_input.ndjson")
BAM_PATH = Path(
    "tests/data/samples_large/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
)
OUTPUT_PATH_NEXTCLADE = OUTPUT_PATH.parent / "nextclade_output.ndjson"


def process(bam_path: Path, output_fp: Path, reference: Path) -> None:
    """Convert BAM to cleartext alignment and write to a ndjson file.

    One Json per read with elements
    - read_id: read identifier
    - query_seq_aligned: the aligned sequence with gaps
    - inserts: list of insertions with elements
            - pos: position of the insertion
            - ins: insertion sequence
    """
    # check if the output directory exists
    output_fp.parent.mkdir(parents=True, exist_ok=True)

    bam_to_cleartext_alignment(bam_path, output_fp, reference)


if __name__ == "__main__":
    N = 1000  # Number of lines to read and process

    process(BAM_PATH, OUTPUT_PATH, REFERENCE_PATH)

    # read the first N lines of the output file
    with open(OUTPUT_PATH, "r") as f:
        lines = [f.readline().replace("'", '"') for _ in range(N)]

    # load the reference sequence
    with open(REFERENCE_PATH, "r") as f:
        reference_seq = f.readlines()
    reference_seq = "".join(reference_seq)

    # Simulate processing times for each sequence
    processing_times = []
    for _ in range(N):
        seq_start_time = time.time()
        # Simulate processing of a sequence
        time.sleep(0.01)  # Simulate a delay for processing
        seq_end_time = time.time()
        processing_times.append(seq_end_time - seq_start_time)

    # Calculate statistics
    mean_time = mean(processing_times)
    variance_time = variance(processing_times)
    min_time = min(processing_times)
    max_time = max(processing_times)

    print(f"Mean processing time: {mean_time:.6f} seconds")
    std_dev_time = variance_time**0.5
    print(f"Standard deviation of processing time: {std_dev_time:.6f} seconds")
    print(f"Minimum processing time: {min_time:.6f} seconds")
    print(f"Maximum processing time: {max_time:.6f} seconds")

    # Estimate time for 5 million reads
    estimated_time_5m_reads = mean_time * 5_000_000
    print(
        f"Estimated time for 5 million reads: {estimated_time_5m_reads / 3600:.2f} hours"
    )

    with open(OUTPUT_PATH_NEXTCLADE, "a") as f:
        for result in lines:
            f.write(json.dumps(result) + "\n")

    print(f"Processed {len(lines)} sequences")
