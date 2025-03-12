"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pysam

from sr2silo.process import paired_end_read_merger
from sr2silo.process.merge import sort_sam_by_qname

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
        paired_end_read_merger(SAM_FP, REF_GENOME_FP, output_merged_sam_fp)

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


def test_sort_sam_by_qname():
    """Test the sort_sam_by_qname function."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        output_sam_path = Path(tmp_dir) / "sorted_output.sam"
        sort_sam_by_qname(SAM_FP, output_sam_path)
        assert output_sam_path.exists(), "Output sam file was not created."
        with pysam.AlignmentFile(str(output_sam_path), "rb") as af:
            so = af.header.get("HD", {}).get("SO")  # type: ignore
            assert so == "queryname", "Output sam file is not sorted by query name."
