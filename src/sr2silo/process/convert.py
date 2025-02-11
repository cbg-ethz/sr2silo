"""Implements conversions from sam to bam and vice versa."""

from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path
from typing import List, Union

import pysam

from sr2silo.process.interface import AAInsertion, Gene

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


def bam_to_fasta(bam_file: Path, fasta_file: Path):
    """
    Convert a BAM file to a FASTA file. Bluntly resolved the sam to fasta.

    Args:
        bam_file: Path to the input BAM file.
        fasta_file: Path to the output FASTQ file.
    """
    # check for proper format
    if not bam_file.suffix.endswith(".bam"):
        raise ValueError("Input file is not a BAM file")
    if not fasta_file.suffix.endswith(".fasta"):
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


# TODO: identify where needed, likely duplicate with sam_to_seq_and_indels
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


def sam_to_seq_and_indels(
    seq: str, cigar: str
) -> Tuple[str, List[AAInsertion], List[Tuple[int, int]]]:
    """
    Processes a SAM file-style sequence and a CIGAR string to return the
    cleartext sequence, along with detailed information about insertions
    and deletions.

    Args:
        seq (str): The sequence string from the SAM file, representing the read.
        cigar (str): The CIGAR string that describes how the sequence aligns
                    to a reference.

    Returns:
        tuple: A tuple containing:
            - cleartext_sequence (str): The sequence aligned to the reference,
                                         excluding insertions and deletions.
            - insertions (list of tuples): A list of tuples, each containing:
                - position (int): The position in the reference where the
                                  insertion occurs.
                - inserted_sequence (str): The sequence that is inserted at
                                           the given position.
            - deletions (list of tuples): A list of tuples, each containing:
                - position (int): The position in the reference where the
                                  deletion starts.
                - length (int): The length of the deletion.

    Example:
        sequence = "AGCTTAGCTAGCTT"
        cigar = "5M1I5M1D3M"
        cleartext, insertions, deletions = sam_to_seq_and_indels(sequence, cigar)

        # Output:
        # Cleartext Sequence: AGCTTAGCTAGC
        # Insertions: [(5, 'A')]
        # Deletions: [(11, 1)]

    Notes:
        - The function assumes that the input sequence and CIGAR string are
          valid and correctly formatted.
        - The CIGAR operations handled include:
            - 'M', '=', 'X': Match or mismatch (aligned to the reference).
            - 'I': Insertion to the reference.
            - 'D': Deletion from the reference.
            - 'N': Skipped region from the reference.
            - 'S': Soft clipping (clipped sequences present in SEQ).
            - 'H': Hard clipping (clipped sequences NOT present in SEQ).
            - 'P': Padding (silent deletion from padded reference).
    """
    parsed_cigar = parse_cigar_new(cigar)
    cleartext_sequence = []
    insertions = []
    deletions = []

    seq_index = 0
    ref_position = 0

    for length, op in parsed_cigar:
        if op == "M" or op == "=" or op == "X":  # Match or mismatch
            cleartext_sequence.append(seq[seq_index : seq_index + length])
            seq_index += length
            ref_position += length
        elif op == "I":  # Insertion to the reference
            insertions.append((ref_position, seq[seq_index : seq_index + length]))
            seq_index += length
        elif op == "D":  # Deletion from the reference
            deletions.append((ref_position, length))
            ref_position += length
        elif op == "N":  # Skipped region from the reference
            ref_position += length
        elif op == "S":  # Soft clipping (clipped sequences present in SEQ)
            seq_index += length
        elif op == "H":  # Hard clipping (clipped sequences NOT present in SEQ)
            pass
        elif op == "P":  # Padding (silent deletion from padded reference)
            pass

    # convert insertions to AAInsertion objects
    insertions = [AAInsertion(position, sequence) for position, sequence in insertions]

    return "".join(cleartext_sequence), insertions, deletions


def get_genes_and_lengths_from_ref(reference_fp: Path) -> Dict[str, Gene]:
    """Load the gene ref fasta and get all the gene names."""
    genes = dict()

    with open(reference_fp, "r") as f:
        awaiting_next_line = False
        for line in f:
            if line.startswith(">"):
                gene = line[1:].strip()
                awaiting_next_line = True
            elif awaiting_next_line:
                reference_length = len(line.strip())
                genes[gene] = Gene(gene, reference_length)
                awaiting_next_line = False
            else:
                continue

    return genes



def sort_and_index_bam(input_bam_fp: Path, output_bam_fp: Path) -> None:
    """
    Function to sort and index the input BAM file,

    implements checks to see if the input BAM file is already sorted and indexed.

    If not sorted and indexed, the function will sort and index the input BAM file.

    """

    # check if sorted and indexed
    if not is_bam_sorted(input_bam_fp) or not is_bam_indexed(
        input_bam_fp
    ):
        logging.info("Sorting and indexing the input BAM file")
        _sort_and_index_bam(
            input_bam_fp, output_bam_fp
        )
    else:
        # copy the input BAM file to the output BAM file
        output_bam_fp.write_bytes(input_bam_fp.read_bytes())
        logging.info("Input BAM file is already sorted and indexed, \
                      copying to output")


def _sort_and_index_bam(input_bam_fp: Path, output_bam_fp: Path) -> None:
    """
    Function to sort and index the input BAM file.

    """

    # Sort the BAM file
    logging.info("Sorting the input BAM file")
    sort_bam_file(input_bam_fp, output_bam_fp)

    # Create index for BAM file if needed
    logging.info("Creating index for the sorted BAM file")
    create_index(output_bam_fp)

    return None


def is_bam_sorted(bam_file):
    """Checks if a BAM file is sorted using pysam.

    Args:
        bam_file (str): Path to the BAM file.

    Returns:
        bool: True if the BAM file is sorted, False otherwise.  Returns None if there's an issue opening the file.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")  # Open in read-binary mode
        is_sorted = bam.header.get("HD", {}).get("SO") == "coordinate"
        bam.close()  # Important: Close the file!
        return is_sorted
    except ValueError as e:
        print(f"Error opening BAM file {bam_file}: {e}")  # Handle file errors
        return None  # Indicate an issue
    except Exception as e:  # Catch other potential errors (like missing index)
        print(f"An unexpected error occurred: {e}")
        return None


def is_bam_indexed(bam_file):
    """Checks if a BAM file has an index (.bai) file.

    Args:bam_file.suffix.endswith
        bam_file (str): Path to the BAM file.

    Returns:
        bool: True if the BAM file has an index, False otherwise. Returns None if there's an issue opening the file.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        has_index = bam.has_index()  # Directly check for index
        bam.close()
        return has_index

    except ValueError as e:
        print(f"Error opening BAM file {bam_file}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None