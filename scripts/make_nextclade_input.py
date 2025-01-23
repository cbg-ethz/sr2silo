"""Make cleartext alignment from BAM file with reference sequence
    for Nextclade input."""

from __future__ import annotations

from pathlib import Path

import nextclade
from sr2silo.process.convert import bam_to_cleartext_alignment

REFERENCE_PATH = Path("resources/sars-cov-2/NC_045512.2.fasta")
OUTPUT_PATH = Path("out/nextclade/nextclade_input.ndjson")
BAM_PATH = Path(
    "tests/data/samples_large/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
)


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
    # process(BAM_PATH, OUTPUT_PATH, REFERENCE_PATH)
    print(nextclade.sum_as_string(1, 3))  # type: ignore
