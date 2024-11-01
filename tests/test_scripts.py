"""This is a sample python file for testing scripts from the source code."""

from __future__ import annotations

import pysam
import subprocess


def bam_to_sam(bam_file):
    """Converts a BAM file to SAM format and returns it as a string.

    Args:
      bam_file: Path to the input BAM file.
    """
    import tempfile

    with tempfile.NamedTemporaryFile(delete=False) as temp_sam:
        with pysam.AlignmentFile(bam_file, "rb") as in_bam, pysam.AlignmentFile(
            temp_sam.name, "w", template=in_bam
        ) as out_sam:
            for read in in_bam:
                out_sam.write(read)
        temp_sam.seek(0)
        return temp_sam.read().decode()


def test_read():
    """Test the read.py script by piping SAM data to it."""
    input_bam = "tests/data/REF_aln_trim_subsample.bam"
    sam_data = bam_to_sam(input_bam)

    # Pipe the SAM data to the script as stdin
    result = subprocess.run(["python", "scripts/dgicev/read.py"], input=sam_data, capture_output=True, text=True)
    assert result.returncode == 0
