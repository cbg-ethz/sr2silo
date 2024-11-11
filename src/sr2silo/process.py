"""This module contains the main functions for processing the data.
"""

import re

from pathlib import Path


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """Parse a cigar string into a list of tuples."""
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")

    parsed_cigar = pattern.findall(cigar)

    return [(op, int(length)) for length, op in parsed_cigar]


def pair_normalize_reads(sam_data: str, output_fasta: Path, output_insertions: Path) -> None:
    """
    Pair and normalize all reads in a SAM file.

    Note that the input SAM file must be read in its entirety before calling this function,
    whilst the output files can be written incrementally.

    Args:
        sam_data: A file-like object containing SAM formatted data.
    Returns:
        A string with merged, normalized reads in FASTA format.
        TODO: annotate the output format
    """
    unpaired = dict()

    with output_fasta.open("w") as fasta_file, output_insertions.open("w") as insertions_file:
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

            if qname in unpaired:
                read1 = unpaired.pop(qname)
                read2 = read

                if read1["pos"] > read2["pos"]:
                    read1, read2 = read2, read1

                index = read1["pos"]
                read1len = len(read1["RESULT_seqUENCE"])
                merged = read1["RESULT_seqUENCE"][: min(read1len, read2["pos"] - read1["pos"])]

                # do deletions cause a problem here?
                gaplen = read1["pos"] + read1len - read2["pos"]
                if gaplen < 0:
                    merged += "N" * (-gaplen)
                    merged += read2["RESULT_seqUENCE"]
                else:
                    overlap_read1 = read1["RESULT_seqUENCE"][read2["pos"] - read1["pos"] :]
                    overlap_read2 = read2["RESULT_seqUENCE"][0 : max(0, gaplen)]

                    overlap_qual1 = read1["RESULT_qual"][read2["pos"] - read1["pos"] :]
                    overlap_qual2 = read2["RESULT_qual"][0 : max(0, gaplen)]

                    # let's set the read1's version by default
                    overlap_result = list(overlap_read1)

                    if overlap_result and overlap_read1 != overlap_read2:
                        if len(overlap_read1) != len(overlap_read2):
                            print("overlaps don't match in size")
                        number_of_diffs = 0
                        for i, (base1, base2) in enumerate(zip(overlap_read1, overlap_read2)):
                            if base1 != base2:
                                # read1 has no quality, and read2 has overlap, so we take it
                                if overlap_qual1[i] == "-" and overlap_read2 != "-":
                                    overlap_result[i] = base2
                                # read2 has better quality, so we take it
                                elif overlap_qual1[i] > overlap_qual2[i]:
                                    overlap_result[i] = base2
                                number_of_diffs += 1

                    merged += "".join(overlap_result) + read2["RESULT_seqUENCE"][max(0, gaplen) :]

                if len(merged) != read2["pos"] + len(read2["RESULT_seqUENCE"]) - read1["pos"]:
                    raise Exception("Length mismatch")

                fasta_file.write(f">{qname}|{read1['pos']}\n{merged}\n")

                merged_insertions = read1["insertions"].copy()
                insertion_index = read1["pos"] + read1len
                merged_insertions += [insert for insert in read2["insertions"] if insert[0] > insertion_index]

                insertions_file.write(f"{qname}\t{merged_insertions}\n")

            else:
                unpaired[qname] = read
        for read_id, unpaired_read in unpaired.items():
            fasta_file.write(f">{read_id}|{unpaired_read['pos']}\n{unpaired_read['RESULT_seqUENCE']}\n")
            insertions_file.write(f"{read_id}\t{unpaired_read['insertions']}\n")
