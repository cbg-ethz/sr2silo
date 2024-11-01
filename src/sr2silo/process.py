"""This module contains the main functions for processing the data.
"""

import re


def parse_cigar(cigar: str) -> list[tuple[str, int]]:
    """Parse a CIGAR string into a list of tuples."""
    pattern = re.compile(r"(\d+)([MIDNSHP=X])")

    parsed_cigar = pattern.findall(cigar)

    return [(op, int(length)) for length, op in parsed_cigar]
