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


def bam_to_cleartext_alignment(
    bam_path: Path, output_fp: Path, reference: Path
) -> None:
    """Convert BAM to cleartext alignment and write to a ndjson file.

    One Json per read with elements
    - read_id: read identifier
    - query_seq_aligned: the aligned sequence with gaps
    - inserts: list of insertions with elements
            - pos: position of the insertion
            - ins: insertion sequence
    """

    # Get the reference length
    logging.info(reference)
    reference_length = 0
    with Path(reference).open() as ref_f:
        for line in ref_f:
            if line.startswith(">"):
                continue
            reference_length += len(line.strip())

    print(f"Reference length: {reference_length}")

    # Ensure the BAM file is indexed
    bam_path_str = str(bam_path)
    sorted_bam_path_str = None
    if not bam_path.with_suffix(".bai").exists():
        try:
            pysam.index(bam_path_str)
        except pysam.SamtoolsError as e:  # type: ignore
            logging.error(f"Error indexing BAM file: {e}")
            sorted_bam_path_str = bam_path_str.replace(".bam", ".sorted.bam")
            pysam.sort("-o", sorted_bam_path_str, bam_path_str)
            pysam.index(sorted_bam_path_str)
            bam_path_str = sorted_bam_path_str

    with Path(output_fp).open("w") as out_f:
        # Open the BAM file
        with pysam.AlignmentFile(bam_path_str, "rb") as samfile:
            for read in samfile.fetch():
                aligned = []  # To store the aligned sequence with gaps
                insertions = []  # To store the insertions
                ref_pos = read.reference_start
                seq_pos = 0

                # validate ciagtuples, query_sequence, and reference_name exist
                if (
                    not read.cigartuples
                    or not read.query_sequence
                    or not read.reference_name
                ):
                    logging.warning(
                        f"Skipping read {read.query_name} due to missing data"
                    )
                    continue

                for operation, length in read.cigartuples:
                    if operation == 0:  # Match (M)
                        aligned.append(read.query_sequence[seq_pos : seq_pos + length])
                        seq_pos += length
                        ref_pos += length
                    elif operation == 1:  # Insertion (I)
                        # Extract insertion sequences
                        insertions.append(
                            (seq_pos, read.query_sequence[seq_pos : seq_pos + length])
                        )
                        aligned.append(
                            read.query_sequence[seq_pos : seq_pos + length]
                        )  # Add insertions to the aligned sequence as well
                        seq_pos += length
                    elif operation == 2:  # Deletion (D)
                        aligned.append(
                            "-" * length
                        )  # Represent deletions as gaps ('-')
                        ref_pos += length
                    elif operation == 4:  # Soft clip (S)
                        seq_pos += length
                    elif operation == 5:  # Hard clip (H)
                        pass  # Don't include hard clipped sequences in the alignment

                # Combine the aligned sequence
                aligned_str = "".join(aligned)

                # Calculate the padding needed for the left and right
                left_padding = "-" * read.reference_start
                right_padding = "-" * (
                    reference_length - len(aligned_str) - read.reference_start
                )

                # Pad the aligned sequence
                padded_alignment = left_padding + aligned_str + right_padding

                # Create a JSON object for the read
                read_json = {
                    "read_id": read.query_name,
                    "query_seq_aligned": padded_alignment,
                    "inserts": (
                        [
                            {"pos": pos + read.reference_start, "ins": list(seq)}
                            for pos, seq in insertions
                        ]
                        if insertions
                        else []
                    ),
                }

                # Write the JSON object to the file
                out_f.write(f"{read_json}\n")

    # Cleanup generated files
    if sorted_bam_path_str:
        Path(sorted_bam_path_str).unlink(missing_ok=True)
        Path(sorted_bam_path_str + ".bai").unlink(missing_ok=True)
    else:
        Path(bam_path_str + ".bai").unlink(missing_ok=True)
