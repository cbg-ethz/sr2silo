"""Implements conversions from sam to bam and vice versa."""

from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path
from typing import List, Union

import pysam

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def sort_bam_file(input_bam_path: Path, output_bam_path: Path):
    """
    Sorts a BAM file using pysam.sort to avoid loading all alignments into memory.

    Args:
        input_bam_path (Path): Path to the input BAM file.
        output_bam_path (Path): Path to the output sorted BAM file.
    """
    try:
        # Convert Path objects to strings for pysam compatibility
        input_bam_str = str(input_bam_path)
        output_bam_str = str(output_bam_path)

        # Using pysam.sort command to sort the BAM file and write to disk incrementally.
        pysam.sort("-o", output_bam_str, input_bam_str)
        print(f"BAM file has been sorted and saved to {output_bam_str}")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise Exception(f"An error occurred: {e}")


def create_index(bam_file: Path):
    """
    Create an index for a BAM file using pysam.

    Args:
        bam_file (str): Path to the input BAM file.
    """
    try:
        # Convert Path object to string for pysam compatibility
        bam_file_str = str(bam_file)

        # Open BamFile and save with 'bai' extension
        pysam.index(bam_file_str)
        print(f"Index created for {bam_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


def bam_to_fasta(bam_file, fasta_file):
    """
    Convert a BAM file to a FASTA file. Bluntly resolved the sam to fasta.

    Args:
        bam_file (str): Path to the input BAM file.
        fasta_file (str): Path to the output FASTQ file.
    """
    # check for proper format
    if not bam_file.endswith(".bam"):
        raise ValueError("Input file is not a BAM file")
    if not fasta_file.endswith(".fasta"):
        raise ValueError("Output file is not a FASTA file")

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(fasta_file, "w") as fq:
            for read in bam.fetch():
                if not read.is_unmapped:
                    name = read.query_name
                    seq = read.query_sequence
                    fq.write(f">{name}\n{seq}\n")


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


def bam_to_fastq_handle_indels(
    bam_file, fastq_file, insertions_file, deletion_char="-"
):
    """
    Convert a BAM file to a FASTQ file, removing insertions and adding a special character for deletions.
    Save the insertions to a separate file. Include alignment positions in the FASTQ file.

    Used to look at the cleartext nucleotide sequence of the reads.

    :param bam_file: Path to the input BAM file
    :param fastq_file: Path to the output FASTQ file
    :param insertions_file: Path to the output file containing insertions
    :param deletion_char: Special character to use for deletions
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(
        fastq_file, "w"
    ) as fastq, open(insertions_file, "w") as insertions:
        for read in bam.fetch():
            if not read.is_unmapped:
                query_sequence = read.query_sequence
                query_qualities = read.query_qualities
                new_sequence = []
                new_qualities = []
                insertion_positions = []

                query_pos = 0
                ref_align_start = read.reference_start
                ref_pos = read.reference_start

                for cigar in read.cigartuples:
                    if cigar[0] == 0:  # Match or mismatch
                        new_sequence.extend(
                            query_sequence[query_pos : query_pos + cigar[1]]
                        )
                        new_qualities.extend(
                            query_qualities[query_pos : query_pos + cigar[1]]
                        )
                        query_pos += cigar[1]
                        ref_pos += cigar[1]
                    elif cigar[0] == 1:  # Insertion
                        insertion_seq = query_sequence[query_pos : query_pos + cigar[1]]
                        insertion_qual = [
                            chr(q + 33)
                            for q in query_qualities[query_pos : query_pos + cigar[1]]
                        ]
                        insertion_positions.append(
                            (ref_pos, insertion_seq, insertion_qual)
                        )
                        query_pos += cigar[1]
                    elif cigar[0] == 2:  # Deletion
                        new_sequence.extend([deletion_char] * cigar[1])
                        new_qualities.extend(
                            [0] * cigar[1]
                        )  # Assigning a low-quality score for deletions
                        ref_pos += cigar[1]

                # Write the modified read to the FASTQ file
                fastq.write(f"@{read.query_name}\n")
                fastq.write(f"{''.join(new_sequence)}\n")
                fastq.write("+\n")
                fastq.write(
                    f"{''.join(chr(q + 33) for q in new_qualities)}\n"
                )  # Phred33 encoding

                # Write the alignment positions to the FASTQ file
                fastq.write(f"alignment_position:{ref_align_start}\n")

                # Write the insertions to the insertions file
                for insertion_pos, insertion_seq, insertion_qual in insertion_positions:
                    insertions.write(
                        f"{read.query_name}\t{insertion_pos}\t{''.join(insertion_seq)}\t{''.join(insertion_qual)}\n"
                    )


def parse_cigar_new(cigar: str) -> List[Tuple[int, str]]:
    """Parse the CIGAR string into a list of tuples."""
    return [
        (int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


# TODO identify where needed and remove this function
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


# TODO: identify where needed, likly duplicate with sam_to_seq_and_indels
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


def pad_alignment(
    sequence: Union[List[str], str],
    reference_start: int,
    reference_length: int,
    unknown_char: str = "N",
) -> str:
    """
    Pad the sequence to match the reference length.

    This function takes a sequence and pads it with a specified character to align it
    with a reference sequence of a given length. The padding is added to both the
    beginning and the end of the sequence as needed.

    Args:
        sequence (Union[List[str], str]): The sequence to be padded.
        reference_start (int): The starting position of the reference sequence.
        reference_length (int): The total length of the reference sequence.
        unknown_char (str, optional): The character to use for padding. Defaults to "N".

    Returns:
        str: The padded sequence as a single string.
    """

    # Combine the aligned sequence
    aligned_str = "".join(sequence)

    # Calculate the padding needed for the left and right
    left_padding = unknown_char * reference_start
    right_padding = unknown_char * (
        reference_length - len(aligned_str) - reference_start
    )

    # Pad the aligned sequence
    padded_alignment = left_padding + aligned_str + right_padding

    return padded_alignment


def bam_to_cleartext_alignment(
    bam_path: Path, output_fp: Path, reference: Path
) -> None:
    """Convert BAM to cleartext alignment and write to a ndjson file.

    Args:
        bam_path: Path to the input BAM file.
        output_fp: Path to the output ndjson.
        reference: Path to the reference sequence as a FASTA file.

    Output:
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

                padded_alignment = pad_alignment(
                    aligned, read.reference_start, reference_length
                )

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
