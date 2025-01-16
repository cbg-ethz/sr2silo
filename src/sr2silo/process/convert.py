"""Implements conversions from sam to bam and vice versa."""

from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path

import pysam


def bam_to_sam(bam_file: Path) -> str:
    """Converts a BAM file to SAM format and returns it as a string.

    Args:
      bam_file: Path to the input BAM file.

    Returns:
      A string containing the SAM format of the input BAM file.
    """

    with tempfile.NamedTemporaryFile(delete=True) as temp_sam:
        with pysam.AlignmentFile(str(bam_file), "rb") as in_bam, pysam.AlignmentFile(
            temp_sam.name, "w", template=in_bam
        ) as out_sam:
            for read in in_bam:
                out_sam.write(read)
        temp_sam.seek(0)
        temp_sam_content = temp_sam.read().decode()
    return temp_sam_content


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """
    Parse a cigar string into a list of tuples.

    Args:
        cigar: A string representing the cigar string.
    Returns:
        A list of tuples where the first element is the operation
        type and the second is the length of the operation.

    Credits: adapted from David Gicev @davidgicev
    """
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")

    parsed_cigar = pattern.findall(cigar)

    return [(op, int(length)) for length, op in parsed_cigar]


def normalize_reads(sam_data: str, output_fasta: Path, output_insertions: Path) -> None:
    """
    Normalize (to clear text sequence using CIGAR)
    all reads in a SAM file output FASTA and insertions files.

    Note that the input SAM file must be read in
    its entirety before calling this function,
    whilst the output files can be written incrementally.

    Args:
        sam_data: A file-like object containing SAM formatted data.
    Returns:
        A string with merged, normalized reads in FASTA format.
        TODO: annotate the output format

    Credits: adapted from David Gicev @davidgicev
    """

    logging.warning(
        "pair_normalize_reads: Nuliotide Insertions are not yet implemented, "
        "{output_insertions} will be empty."
    )

    unpaired = dict()

    with output_fasta.open("w") as fasta_file, output_insertions.open(
        "w"
    ) as insertions_file:
        for line in sam_data.splitlines():
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")

            qname = fields[0]  # Query template NAME
            pos = int(fields[3])  # 1-based leftmost mapping position
            cigar = parse_cigar(fields[5])  # cigar string
            seq = fields[9]  # segment sequence
            qual = fields[10]  # ASCII of Phred-scaled base quality + 33

            result_sequence = ""
            result_qual = ""
            index = 0
            inserts = []

            for operation in cigar:
                ops_type, count = operation
                if ops_type == "S":
                    index += count
                    continue
                if ops_type == "M":
                    result_sequence += seq[index : index + count]
                    result_qual += qual[index : index + count]
                    index += count
                    continue
                if ops_type == "D":
                    result_sequence += "-" * count
                    result_qual += "!" * count
                    continue
                if ops_type == "I":
                    inserts.append((index + pos, seq[index : index + count]))
                    index += count
                    continue

            read = {
                "pos": pos,
                "cigar": cigar,
                "RESULT_seqUENCE": result_sequence,
                "RESULT_qual": result_qual,
                "insertions": inserts,
            }

            fasta_file.write(f">{qname}|{read['pos']}\n{read['RESULT_seqUENCE']}\n")

            insertions = read["insertions"].copy()
            # insertion_index = read["pos"] + len(read["RESULT_seqUENCE"])

            insertions_file.write(f"{qname}\t{insertions}\n")

        for read_id, unpaired_read in unpaired.items():
            fasta_file.write(
                f">{read_id}|{unpaired_read['pos']}\n{unpaired_read['RESULT_seqUENCE']}\n"
            )
            insertions_file.write(f"{read_id}\t{unpaired_read['insertions']}\n")
