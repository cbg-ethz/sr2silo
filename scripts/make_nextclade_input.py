"""Make cleartext alignment from BAM file with reference sequence
    for Nextclade input."""

from __future__ import annotations

import json
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
    process(BAM_PATH, OUTPUT_PATH, REFERENCE_PATH)
    # read in the first line of the output file
    with open(OUTPUT_PATH, "r") as f:
        first_line = f.readline()
    # parsel as first line of ndjson file, which is a dictionary
    # handle : json.decoder.JSONDecodeError: Expecting property name enclosed in double quotes: line 1 column 2 (char 1)
    first_line = first_line.replace("'", '"')
    first_line_dict = json.loads(first_line)
    print(first_line_dict.keys())

    # load the reference sequence
    with open(REFERENCE_PATH, "r") as f:
        reference_seq = f.readlines()
    reference_seq = "".join(reference_seq[1:]).strip()
    # make one line
    reference_seq = reference_seq.replace("\n", "")

    # test nextclade
    qry_seq = first_line_dict["query_seq_aligned"]
    ref_seq = reference_seq

    # check that both sequences are the same length
    print(
        f"Reference and Aligned Query Seq have the same length {len(qry_seq) == len(ref_seq)}"
    )

    # gene ref path
    GENE_MAP_GFF = "nextclade/data/sars-cov-2/genemap.gff"

    # print(reference_seq)
    print(nextclade.translate_aa_align(ref_seq, qry_seq, GENE_MAP_GFF))  # type: ignore
