"""Implements conversions from sam to bam and vice versa."""

import tempfile
import pysam


def bam_to_sam(bam_file):
    """Converts a BAM file to SAM format and returns it as a string.

    Args:
      bam_file: Path to the input BAM file.
    """

    with tempfile.NamedTemporaryFile(delete=True) as temp_sam:
        with pysam.AlignmentFile(bam_file, "rb") as in_bam, pysam.AlignmentFile(
            temp_sam.name, "w", template=in_bam
        ) as out_sam:
            for read in in_bam:
                out_sam.write(read)
        temp_sam.seek(0)
        temp_sam_content = temp_sam.read().decode()
    return temp_sam_content
