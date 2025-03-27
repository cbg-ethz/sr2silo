"""Implements conversions from sam to bam and vice versa."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Tuple, Union

import pysam

from sr2silo.process.interface import (
    Gene,
    GeneName,
    GeneSet,
    Insertion,
)

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def get_sequence_from_fasta(fasta_fp: Path) -> str:
    """
    Extracts sequence from a FASTA file, ignoring headers.

    Args:
        fasta_fp (Path): Path to the FASTA file.

    Returns:
        str: The concatenated sequence without headers.

    Notes:
        If multiple sequences are in the file, they will be concatenated.
    """
    sequence = ""
    with open(fasta_fp, "r") as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and header lines (starting with '>')
            if not line or line.startswith(">"):
                continue
            sequence += line
    return sequence


def sort_bam_file(
    input_bam_path: Path, output_bam_path: Path, sort_by_qname: bool = False
):
    """
    Sorts a BAM file using pysam.sort by alignment positions by default,
    but can also sort by query name if specified.

    Args:
        input_bam_path (Path): Path to the input BAM file.
        output_bam_path (Path): Path to the output sorted BAM file.
        sort_by_qname (bool, optional): If True, sorts by query name. Defaults to False.
    """
    try:
        # Convert Path objects to strings for pysam compatibility
        input_bam_str = str(input_bam_path)
        output_bam_str = str(output_bam_path)

        # Build sort arguments based on sorting option.
        if sort_by_qname:
            # Using the -n flag to sort by query name.
            pysam.sort("-n", "-o", output_bam_str, input_bam_str)
            logging.info(
                f"BAM file has been sorted by query name and saved to {output_bam_str}"
            )
        else:
            pysam.sort("-o", output_bam_str, input_bam_str)
            logging.info(
                f"BAM file has been sorted by coordinate and saved to {output_bam_str}"
            )
    except Exception as e:
        print(f"An error occurred: {e}")
        raise Exception(f"An error occurred: {e}")


def create_index(bam_file: Path):
    """
    Create an index for a BAM file using pysam.

    Leaves a .bai file in the same directory as the input BAM file.

    Args:
        bam_file (str): Path to the input BAM file.
    """
    try:
        # Convert Path object to string for pysam compatibility
        bam_file_str = str(bam_file)

        # Open BamFile and save with 'bai' extension
        pysam.index(bam_file_str)
        logging.info(f"Index created for {bam_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


def bam_to_fasta_query(bam_file: Path, fasta_file: Path):
    """
    Convert a BAM file to a FASTA file. Bluntly resolved the sam to fasta,
    removing soft clippings, keeping insertions and deletions, ignoring skipped
    regions and paddings.

    Outputs the sequence as it is read from the molecule.

    Args:
        bam_file: Path to the input BAM file.
        fasta_file: Path to the output FASTQ file.
    """
    # check for proper format
    if not bam_file.suffix.endswith(".bam"):
        raise ValueError("Input file is not a BAM file")
    if not fasta_file.suffix.endswith(".fasta"):
        raise ValueError("Output file is not a FASTA file")

    # check if index exists, make index if not
    if not bam_file.with_suffix(".bai").exists():
        sort_and_index_bam(bam_file, bam_file)

    with pysam.AlignmentFile(str(bam_file), "rb") as bam:
        with open(fasta_file, "w") as fq:
            for read in bam.fetch():
                if not read.is_unmapped:
                    name = read.query_name
                    full_seq = read.query_sequence

                    # Remove soft clipping
                    if full_seq is None:
                        continue

                    # Get CIGAR string to identify soft clippings
                    cigar_tuples = read.cigartuples
                    if cigar_tuples is None:
                        continue

                    # Calculate sequence without soft clippings
                    start = 0
                    end = len(full_seq)

                    # Check for soft clipping at start (CIGAR op 4 = S)
                    if cigar_tuples[0][0] == 4:
                        start = cigar_tuples[0][1]

                    # Check for soft clipping at end
                    if cigar_tuples[-1][0] == 4:
                        end -= cigar_tuples[-1][1]

                    # Extract sequence without soft clippings
                    seq = full_seq[start:end]

                    fq.write(f">{name}\n{seq}\n")


def bam_to_sam(bam_file: Path, sam_file: Path) -> None:
    """Converts a BAM file to SAM format and writes it to the specified output file.

    Args:
      bam_file: Path to the input BAM file.
      sam_file: Path to the output SAM file.
    """
    with pysam.AlignmentFile(str(bam_file), "rb") as in_bam, pysam.AlignmentFile(
        str(sam_file), "w", template=in_bam
    ) as out_sam:
        for read in in_bam:
            out_sam.write(read)
    logging.info(f"BAM file {bam_file} has been converted to SAM file {sam_file}")


def sam_to_bam(sam_file: Path, bam_file: Path):
    """
    Convert a SAM file to a BAM file.

    Args:
        sam_file (Path): Path to the input SAM file.
        bam_file (Path): Path to the output BAM file.
    """
    # check for proper format
    if not sam_file.suffix.endswith(".sam"):
        raise ValueError("Input file is not a SAM file")
    if not bam_file.suffix.endswith(".bam"):
        raise ValueError("Output file is not a BAM file")

    with pysam.AlignmentFile(str(sam_file), "r") as sam:
        with pysam.AlignmentFile(str(bam_file), "wb", header=sam.header) as bam:
            for read in sam:
                bam.write(read)

    logging.info(f"SAM file {sam_file} has been converted to BAM file {bam_file}")


def bam_to_fastq_handle_indels(
    bam_file: Path,
    out_fastq_fp: Path,
    out_insertions_fp: Path,
    deletion_char: str = "-",
    skipped_char: str = "N",
):
    """
    Convert a BAM file to a FASTQ file, removing insertions and adding a
    special character for deletions, skipped regions, and soft clipping.
    Save the insertions to a separate file.
    Include alignment positions in the FASTQ file.

    Coordinates are 1-based.

    :param bam_file: Path to the input BAM file
    :param out_fastq_fp: Path to the output FASTQ file
    :param out_insertions_fp: Path to the output file containing insertions
    :param deletion_char: Special character to use for deletions/skipped regions
    :param skipped_char: Special character to use for skipped regions
    """
    with pysam.AlignmentFile(str(bam_file), "rb") as bam, open(
        out_fastq_fp, "w"
    ) as fastq, open(out_insertions_fp, "w") as insertions:
        for read in bam.fetch():
            if not read.is_unmapped:
                logging.debug(f"Processing read: {read.query_name}")
                query_sequence = read.query_sequence if read.query_sequence else ""
                query_qualities = read.query_qualities if read.query_qualities else ""
                new_sequence = []
                new_qualities = []
                insertion_positions = []

                query_pos = 0
                ref_align_start = read.reference_start
                ref_pos = read.reference_start

                if read.cigartuples is None:
                    continue

                for cigar in read.cigartuples:

                    # Handle the CIGAR operations
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
                            chr(int(q) + 33)
                            for q in query_qualities[query_pos : query_pos + cigar[1]]
                        ]
                        insertion_positions.append(
                            (ref_pos, insertion_seq, insertion_qual)
                        )
                        query_pos += cigar[1]
                    elif cigar[0] == 2:  # Deletion
                        new_sequence.extend([deletion_char] * cigar[1])
                        new_qualities.extend([0] * cigar[1])
                        ref_pos += cigar[1]
                    elif cigar[0] == 3:  # Skipped region from the reference
                        new_sequence.extend([skipped_char] * cigar[1])
                        new_qualities.extend([0] * cigar[1])
                        ref_pos += cigar[1]
                    elif cigar[0] == 4:  # Soft clipping
                        # Skip soft clipped bases (they are not aligned)
                        query_pos += cigar[1]
                    elif cigar[0] == 5:  # Hard clipping
                        # Hard clipped bases are not present in the read sequence,
                        # so no update to query_pos is required.
                        pass
                    elif cigar[0] == 6:  # Padding
                        # Padding is a silent deletion from padded reference
                        # No action needed for this case.
                        pass

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


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse the CIGAR string into a list of tuples."""
    return [
        (int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


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
) -> Tuple[str, List[Insertion], List[Tuple[int, int]]]:
    """
    Processes a SAM file-style sequence (nuclitide / amino acids) and a CIGAR
    string to return the cleartext sequence, along with detailed information
    about insertions and deletions.

    Args:
        seq (str): The sequence string from the SAM file, representing the read.
        cigar (str): The CIGAR string that describes how the sequence aligns
                    to a reference.

    Returns:
        tuple: A tuple containing:
            - cleartext_sequence (str): The sequence aligned to the reference,
                                         excluding insertions and deletions.
            - insertions (list of Insertion): A list of Insertion objects
            - deletions (list of tuples): A list of tuples, each containing:
                - position (int): The position in the reference where the
                                  deletion starts.
                - length (int): The length of the deletion.

    Example:
        sequence = "AGCTTAGCTAGCTT"
        cigar = "5M1I5M1D3M"

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
    parsed_cigar = parse_cigar(cigar)
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

    insertions = [
        Insertion(position=ins_pos, sequence=ins_seq) for ins_pos, ins_seq in insertions
    ]

    return "".join(cleartext_sequence), insertions, deletions


def get_gene_set_from_ref(reference_fp: Path) -> GeneSet:
    """Load the gene ref fasta and create a GeneSet with gene short
    names and lengths."""
    genes = dict()
    with open(reference_fp, "r") as f:
        # Strip lines and ignore blank lines
        lines = [line.strip() for line in f if line.strip()]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            gene_name = GeneName(lines[i][1:].strip())
            sequence = ""
            # If the next line exists and is not a header, use it as the sequence
            if i + 1 < len(lines) and not lines[i + 1].startswith(">"):
                sequence = lines[i + 1].strip()
                i += 2
            else:
                i += 1
            if sequence:
                genes[gene_name] = Gene(gene_name, len(sequence))
        else:
            i += 1
    gene_set = GeneSet(list(genes.values()))
    # Optional warnings
    if not gene_set.get_gene_name_list():
        logging.warning("No genes found in the reference file")
    if len(gene_set.get_gene_name_list()) < len(lines) / 2:
        logging.warning("Some genes were skipped in the reference file")
        logging.warning("Parsed genes: %s", gene_set.get_gene_name_list())
    return gene_set


def sort_and_index_bam(input_bam_fp: Path, output_bam_fp: Path) -> None:
    """
    Function to sort and index the input BAM file,

    implements checks to see if the input BAM file is already sorted and indexed.

    If not sorted and indexed, the function will sort and index the input BAM file.

    """

    # check if sorted and indexed
    if not is_bam_sorted(input_bam_fp) or not is_bam_indexed(input_bam_fp):
        logging.info("Sorting and indexing the input BAM file")
        _sort_and_index_bam(input_bam_fp, output_bam_fp)
    else:
        # copy the input BAM file to the output, along with the index
        output_bam_fp.write_bytes(input_bam_fp.read_bytes())
        # note then ending is .bam.bai
        output_bam_fp.with_suffix(".bam.bai").write_bytes(
            input_bam_fp.with_suffix(".bam.bai").read_bytes()
        )
        logging.info(
            "Input BAM file is already sorted and indexed, \
                      copying to output"
        )


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
    """Checks if a BAM file is sorted by genomic coordinates using pysam.

    Args:
        bam_file (str): Path to the BAM file.

    Returns:
        bool: True if the BAM file is sorted, False otherwise.
        Returns None if there's an issue opening the file.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")  # Open in read-binary mode
        is_sorted = (
            bam.header.get("HD", {}).get("SO") == "coordinate"  # pyright: ignore
        )  # pyright: ignore
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

    Args:
        bam_file (str): Path to the BAM file.

    Returns:
        bool: True if the BAM file has an index, False otherwise.
        Returns None if there's an issue opening the file.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        has_index = bam.has_index()
        bam.close()
        return has_index

    except ValueError as e:
        print(f"Error opening BAM file {bam_file}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None


def split_bam(
    input_bam: Path, out_dir: Path, chunk_size: int, prefix="split_"
) -> List[Path]:
    """
    Split a BAM file into smaller BAM files, each containing up to chunk_size reads.

    Parameters:
        input_bam (str): Path to the input BAM file.
        chunk_size (int): Number of reads per output BAM file.
        prefix (str): Prefix for the output BAM files (default: "split_").

    Returns:
        list: List of paths to the output BAM files.
    """
    # Open the input BAM file for reading in binary mode
    bamfile = pysam.AlignmentFile(str(input_bam), "rb")

    chunk_num = 1  # Initialize chunk number
    count = 0  # Initialize read counter
    current_chunk_file = None  # Initialize current chunk file as None

    # Iterate through each read in the input BAM file
    for read in bamfile:
        # If count is 0, start a new chunk file
        if count == 0:
            # Close the previous chunk file if it exists
            if current_chunk_file is not None:
                current_chunk_file.close()
            # Open a new chunk file with the original header
            current_chunk_file = pysam.AlignmentFile(
                str(out_dir / f"{prefix}{chunk_num}.bam"), "wb", header=bamfile.header
            )
            chunk_num += 1  # Increment chunk number for the next file

        # Write the current read to the chunk file
        if current_chunk_file is not None:
            current_chunk_file.write(read)
            count += 1  # Increment read counter

        # If the chunk size is reached, reset the counter
        if count == chunk_size:
            count = 0

    # Close the last chunk file if it was opened
    if current_chunk_file is not None:
        current_chunk_file.close()

    # Close the input BAM file
    bamfile.close()

    return list(out_dir.glob(f"{prefix}*.bam"))


def is_sorted_qname(bam_file):
    """Checks if a BAM/SAM file is sorted by QNAME using pysam.

    Args:
        bam_file (str): Path to the BAM/SAM file.

    Returns:
        bool: True if the file is sorted by query name, False otherwise.
        None if there's an issue opening the file.
    """
    # Choose file mode based on extension (assumes .sam files are uncompressed)
    mode = "rb" if Path(bam_file).suffix == ".bam" else "r"
    try:
        with pysam.AlignmentFile(bam_file, mode) as af:
            # Header HD should define SO as "queryname" for a QNAME sorted file.
            so = af.header.get("HD", {}).get("SO")  # type: ignore
            return so == "queryname"
    except ValueError as e:
        print(f"Error opening file {bam_file}: {e}")
        return None
    except Exception as e:  # pragma: no cover
        print(f"An unexpected error occurred: {e}")
        return None


def had_SQ_header(sam_file):
    """Checks if a SAM/BAM file has @SQ headers using pysam.

    Args:
        sam_file (str): Path to the SAM/BAM file.

    Returns:
        bool: True if the file contains @SQ headers, False otherwise.
    """
    try:
        with pysam.AlignmentFile(sam_file, "r") as af:
            # Check if there is an 'SQ' entry in the header
            return "SQ" in af.header and len(af.header["SQ"]) > 0  # type: ignore
    except Exception as e:  # pragma: no cover
        print(f"Error reading SAM file {sam_file}: {e}")
        return False
