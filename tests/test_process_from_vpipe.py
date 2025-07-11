"""Tests for the import_to_loculus module, specifically for the skip_merge option."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest

from sr2silo.process_from_vpipe import nuc_align_to_silo_njson


@pytest.mark.parametrize("skip_merge", [True, False])
def test_skip_merge_option(
    real_sample_files_import_to_loculus: dict[str, Any],
    tmp_path: Path,
    skip_merge: bool,
):
    """Test that the skip_merge option correctly bypasses the read pair merging step.

    Args:
        real_sample_files_import_to_loculus: Fixture providing sample input files
        tmp_path: Pytest fixture for temporary directory
        skip_merge: Whether to skip the read pair merging step
    """
    sample_files = real_sample_files_import_to_loculus
    output_fp = tmp_path / "output.ndjson.zst"

    # Create a test mock directory and file to avoid file not found errors
    mock_result_dir = tmp_path / "results"
    mock_result_dir.mkdir(parents=True, exist_ok=True)
    mock_merged_sam = mock_result_dir / "mock_merged.sam"
    mock_merged_sam.write_text("# Mock SAM content")

    # We'll patch multiple functions to isolate the test and prevent file system errors
    with patch(
        "sr2silo.process_from_vpipe.paired_end_read_merger"
    ) as mock_merger, patch(
        "sr2silo.process_from_vpipe.parse_translate_align_in_batches",
        return_value=output_fp,
    ) as _mock_parse, patch(  # ruff: noqa: E501
        "sr2silo.process_from_vpipe.Path.mkdir", return_value=None
    ) as _mock_mkdir, patch(
        "sr2silo.process_from_vpipe.sort_bam_file"
    ) as _mock_sort, patch(
        "sr2silo.process_from_vpipe.bam_to_sam"
    ) as _mock_bam_to_sam, patch(
        "sr2silo.process_from_vpipe.sam_to_bam"
    ) as _mock_sam_to_bam, patch(
        "sr2silo.process_from_vpipe.Path.unlink"
    ) as _mock_unlink, patch(
        "pathlib.Path.exists", return_value=True
    ), patch(
        "sr2silo.process_from_vpipe.Path.parent", new_callable=MagicMock
    ) as mock_parent:

        # Set up the mock for Path.parent to return a Path that contains a
        # 'results' directory
        mock_parent_path = MagicMock()
        mock_parent_path.__truediv__.return_value = mock_result_dir
        mock_parent.return_value = mock_parent_path

        # Call the function with skip_merge=True or False
        nuc_align_to_silo_njson(
            input_file=sample_files["input_file"],
            sample_id=sample_files["sample_id"],
            timeline_file=sample_files["timeline_file"],
            output_fp=output_fp,
            nuc_ref_fp=sample_files["nuc_ref_fp"],
            aa_ref_fp=sample_files["aa_ref_fp"],
            skip_merge=skip_merge,
        )

        # Check if paired_end_read_merger was called based on skip_merge value
        if skip_merge:
            # If skip_merge is True, paired_end_read_merger should NOT be called
            mock_merger.assert_not_called()
        else:
            # If skip_merge is False, paired_end_read_merger SHOULD be called
            assert mock_merger.call_count == 1


def test_skip_merge_file_handling(
    real_sample_files_import_to_loculus: dict[str, Any], tmp_path: Path
):
    """Test file handling when skip_merge is True.

    Verify that when skip_merge is True:
    1. For BAM files, the input file is used directly
    2. For SAM files, the file is converted to BAM format
    """
    sample_files = real_sample_files_import_to_loculus
    output_fp = tmp_path / "output.ndjson.zst"

    # Create test files and mocks
    mock_result_dir = tmp_path / "results"
    mock_result_dir.mkdir(parents=True, exist_ok=True)

    # Create a test SAM file to test the SAM->BAM conversion case
    test_sam_path = tmp_path / "test.sam"
    test_sam_path.write_text("# Test SAM file")

    # Set up common mocks
    with patch(
        "sr2silo.process_from_vpipe.parse_translate_align_in_batches",
        return_value=output_fp,
    ), patch("sr2silo.process_from_vpipe.Path.mkdir", return_value=None), patch(
        "pathlib.Path.exists", return_value=True
    ), patch(
        "sr2silo.process_from_vpipe.Path.parent", new_callable=MagicMock
    ) as mock_parent:

        # Set up the mock for Path.parent to return a Path that contains a
        # 'results' directory
        mock_parent_path = MagicMock()
        mock_parent_path.__truediv__.return_value = mock_result_dir
        mock_parent.return_value = mock_parent_path

        # Test with a BAM file
        with patch("sr2silo.process_from_vpipe.sam_to_bam") as mock_sam_to_bam:
            nuc_align_to_silo_njson(
                input_file=sample_files["input_file"],  # This is a BAM file
                sample_id=sample_files["sample_id"],
                timeline_file=sample_files["timeline_file"],
                output_fp=output_fp,
                nuc_ref_fp=sample_files["nuc_ref_fp"],
                aa_ref_fp=sample_files["aa_ref_fp"],
                skip_merge=True,
            )

            # For BAM files with skip_merge=True, sam_to_bam should not be called
            # to convert the input file
            assert mock_sam_to_bam.call_count == 0

        # Test with a SAM file
        with patch("sr2silo.process_from_vpipe.sam_to_bam") as mock_sam_to_bam:
            nuc_align_to_silo_njson(
                input_file=test_sam_path,  # This is a SAM file
                sample_id=sample_files["sample_id"],
                timeline_file=sample_files["timeline_file"],
                output_fp=output_fp,
                nuc_ref_fp=sample_files["nuc_ref_fp"],
                aa_ref_fp=sample_files["aa_ref_fp"],
                skip_merge=True,
            )

            # For SAM files with skip_merge=True, sam_to_bam should be called once
            # to convert the input file to BAM format
            assert mock_sam_to_bam.call_count == 1
