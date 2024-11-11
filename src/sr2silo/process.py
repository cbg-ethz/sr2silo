"""This module contains the main functions for processing the data.
"""

import re

from pathlib import Path


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """Parse a CIGAR string into a list of tuples."""
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

            QNAME = fields[0]  # Query template NAME
            POS = int(fields[3])  # 1-based leftmost mapping POSition
            CIGAR = parse_cigar(fields[5])  # CIGAR string
            SEQ = fields[9]  # segment SEQuence
            QUAL = fields[10]  # ASCII of Phred-scaled base QUALity + 33

            result_sequence = ""
            result_qual = ""
            index = 0
            inserts = []

            for operation in CIGAR:
                type, count = operation
                if type == "S":
                    index += count
                    continue
                if type == "M":
                    result_sequence += SEQ[index : index + count]
                    result_qual += QUAL[index : index + count]
                    index += count
                    continue
                if type == "D":
                    result_sequence += "-" * count
                    result_qual += "!" * count
                    continue
                if type == "I":
                    inserts.append((index + POS, SEQ[index : index + count]))
                    index += count
                    continue

            read = {
                "POS": POS,
                "CIGAR": CIGAR,
                "RESULT_SEQUENCE": result_sequence,
                "RESULT_QUAL": result_qual,
                "insertions": inserts,
            }

            if QNAME in unpaired:
                read1 = unpaired.pop(QNAME)
                read2 = read

                if read1["POS"] > read2["POS"]:
                    read1, read2 = read2, read1

                index = read1["POS"]
                read1len = len(read1["RESULT_SEQUENCE"])
                merged = read1["RESULT_SEQUENCE"][: min(read1len, read2["POS"] - read1["POS"])]

                # do deletions cause a problem here?
                gaplen = read1["POS"] + read1len - read2["POS"]
                if gaplen < 0:
                    merged += "N" * (-gaplen)
                    merged += read2["RESULT_SEQUENCE"]
                else:
                    overlap_read1 = read1["RESULT_SEQUENCE"][read2["POS"] - read1["POS"] :]
                    overlap_read2 = read2["RESULT_SEQUENCE"][0 : max(0, gaplen)]

                    overlap_qual1 = read1["RESULT_QUAL"][read2["POS"] - read1["POS"] :]
                    overlap_qual2 = read2["RESULT_QUAL"][0 : max(0, gaplen)]

                    # let's set the read1's version by default
                    overlap_result = list(overlap_read1)

                    if len(overlap_result) and overlap_read1 != overlap_read2:
                        # print("", QNAME)
                        if len(overlap_read1) != len(overlap_read2):
                            print("overlaps don't match in size")
                        number_of_diffs = 0
                        for i in range(len(overlap_read1)):
                            if overlap_read1[i] != overlap_read2[i]:
                                if overlap_qual1[i] == "-" and overlap_read2 != "-":
                                    overlap_result[i] = overlap_read2[i]
                                if overlap_qual1[i] > overlap_qual2[i]:
                                    overlap_result[i] = overlap_read2[i]
                                # print("diff in position ", i, ": ", overlap_read1[i], "/", overlap_read2[i])
                                number_of_diffs += 1

                    merged += "".join(overlap_result) + read2["RESULT_SEQUENCE"][max(0, gaplen) :]

                if len(merged) != read2["POS"] + len(read2["RESULT_SEQUENCE"]) - read1["POS"]:
                    raise Exception("Length mismatch")

                fasta_file.write(f">{QNAME}|{read1['POS']}\n{merged}\n")

                merged_insertions = read1["insertions"].copy()
                insertion_index = read1["POS"] + read1len
                merged_insertions += [insert for insert in read2["insertions"] if insert[0] > insertion_index]

                insertions_file.write(f"{QNAME}\t{merged_insertions}\n")

            else:
                unpaired[QNAME] = read
        for id in unpaired:
            fasta_file.write(f">{id}|{unpaired[id]['POS']}\n{unpaired[id]['RESULT_SEQUENCE']}\n")
            insertions_file.write(f"{id}\t{unpaired[id]['insertions']}\n")
