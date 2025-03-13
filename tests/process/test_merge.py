"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import tempfile
from pathlib import Path

from sr2silo.process import paired_end_read_merger

from sr2silo.process.convert import sort_sam_by_qname

SAM_FP = Path("tests/data/REF_aln_trim_subsample_expected.sam")
REF_GENOME_FP = Path("resources/sars-cov-2/nuc_reference_genomes.fasta")


def test_paired_end_read_merger():
    """Test the paired_end_read_merger function.

    test for:
    - read ids are unique
    - execute the function without errors
    """

    with tempfile.TemporaryDirectory() as tmp_dir:
        output_merged_sam_fp = Path(tmp_dir) / "merged.sam"

        # Test with unsorted SAM file (should fail)
        try:
            paired_end_read_merger(SAM_FP, REF_GENOME_FP, output_merged_sam_fp)
        except ValueError as e:
            assert "not sorted by QNAME" in str(e)

        # Sort the SAM file
        sorted_sam_fp = Path(tmp_dir) / "sorted.sam"
        sort_sam_by_qname(SAM_FP, sorted_sam_fp)

        # Test with sorted SAM file (should pass)
        paired_end_read_merger(sorted_sam_fp, REF_GENOME_FP, output_merged_sam_fp)

        with open(output_merged_sam_fp, "r") as f:
            lines = f.readlines()

        # assert that each read id is unique
        read_ids = set()
        for line in lines:
            if line.startswith("@"):
                continue
            read_id = line.split("\t")[0]
            assert read_id not in read_ids
            read_ids.add(read_id)


def test_paired_end_read_merger_exceptions():
    """Test the exceptions in paired_end_read_merger function."""

    with tempfile.TemporaryDirectory() as tmp_dir:
        output_merged_sam_fp = Path(tmp_dir) / "merged.sam"

        # Test FileNotFoundError for nuc_align_sam_fp
        try:
            paired_end_read_merger(
                Path("non_existent.sam"), REF_GENOME_FP, output_merged_sam_fp
            )
        except FileNotFoundError as e:
            assert "File not found" in str(e)

        # Test FileNotFoundError for ref_genome_fasta_fp
        try:
            paired_end_read_merger(
                SAM_FP, Path("non_existent.fasta"), output_merged_sam_fp
            )
        except FileNotFoundError as e:
            assert "File not found" in str(e)

        # Test ValueError for missing @SQ headers
        with tempfile.NamedTemporaryFile(delete=False) as tmp_sam:
            tmp_sam.write(b"@HD\tVN:1.0\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n")
            tmp_sam_fp = Path(tmp_sam.name)

        # Sort the temporary SAM file
        sorted_tmp_sam_fp = Path(tmp_dir) / "sorted_tmp.sam"
        sort_sam_by_qname(tmp_sam_fp, sorted_tmp_sam_fp)

        try:
            paired_end_read_merger(
                sorted_tmp_sam_fp, REF_GENOME_FP, output_merged_sam_fp
            )
        except ValueError as e:
            assert "does not have @SQ headers" in str(e)

        # Clean up temporary file
        tmp_sam_fp.unlink()
