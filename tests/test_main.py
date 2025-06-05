"""Tests for the main CLI application of sr2silo.
These tests cover the command-line interface and its
functionality, ensuring that the application behaves
as expected when invoked with various commands and options.
"""

from __future__ import annotations

import os

from typer.testing import CliRunner

from sr2silo.main import app

runner = CliRunner()


def test_help():
    """Test that help command works and shows correct information."""
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "sr2silo" in result.stdout
    assert "Convert Short-Read nucleotide .bam alignments" in result.stdout


def test_no_args():
    """Test behavior when no arguments are provided."""
    result = runner.invoke(app)
    assert result.exit_code == 0
    assert "Well, you gotta decide what to do" in result.stdout


def test_run_command():
    """Test the run subcommand."""
    result = runner.invoke(app, ["run"])
    assert result.exit_code == 0
    assert "Not yet implemented" in result.stdout


def test_process_from_vpipe_help():
    """Test the help output for process-from-vpipe command."""
    result = runner.invoke(app, ["process-from-vpipe", "--help"])
    assert result.exit_code == 0
    assert "V-PIPE to SILO conversion" in result.stdout


def test_process_from_vpipe_missing_required():
    """Test process-from-vpipe fails when required arguments are missing."""
    result = runner.invoke(app, ["process-from-vpipe"])
    assert result.exit_code != 0  # Should fail due to missing required options


def test_process_from_vpipe_with_real_files(real_sample_files_import_to_loculus):
    """Test process-from-vpipe with real sample files."""

    result = runner.invoke(
        app,
        [
            "process-from-vpipe",
            "--input-file",
            str(real_sample_files_import_to_loculus["input_file"]),
            "--sample-id",
            real_sample_files_import_to_loculus["sample_id"],
            "--batch-id",
            real_sample_files_import_to_loculus["batch_id"],
            "--timeline-file",
            str(real_sample_files_import_to_loculus["timeline_file"]),
            "--primer-file",
            str(real_sample_files_import_to_loculus["primer_file"]),
            "--output-fp",
            str(real_sample_files_import_to_loculus["output_file"]),
            "--reference",
            real_sample_files_import_to_loculus["reference"],
        ],
    )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout


def test_submit_to_loculus_command_help():
    """Test the help output for submit-to-loculus command."""
    result = runner.invoke(app, ["submit-to-loculus", "--help"])
    assert result.exit_code == 0
    assert "Upload processed file to S3 and submit to SILO/Loculus" in result.stdout


def test_submit_to_loculus_command_missing_required():
    """Test submit-to-loculus command fails when required arguments are missing."""
    result = runner.invoke(app, ["submit-to-loculus"])
    assert result.exit_code != 0  # Should fail due to missing required options


def test_submit_to_loculus_command_nonexistent_file():
    """Test submit-to-loculus command fails when processed file doesn't exist."""
    result = runner.invoke(
        app,
        [
            "submit-to-loculus",
            "--processed-file",
            "/tmp/nonexistent.ndjson.zst",
            "--sample-id",
            "test_sample",
        ],
    )
    assert result.exit_code == 1


def test_submit_to_loculus_command_wrong_extension():
    """Test submit-to-loculus command fails when file has wrong extension."""
    # Create a temporary file with wrong extension
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as tmp:
        tmp.write(b"test content")
        tmp_path = tmp.name

    try:
        result = runner.invoke(
            app,
            [
                "submit-to-loculus",
                "--processed-file",
                tmp_path,
                "--sample-id",
                "test_sample",
            ],
        )
        assert result.exit_code == 1
    finally:
        import os

        os.unlink(tmp_path)


def test_process_from_vpipe_with_skip_merge(real_sample_files_import_to_loculus):
    """Test process-from-vpipe with skip-merge option to bypass read pair merging."""
    result = runner.invoke(
        app,
        [
            "process-from-vpipe",
            "--input-file",
            str(real_sample_files_import_to_loculus["input_file"]),
            "--sample-id",
            real_sample_files_import_to_loculus["sample_id"],
            "--batch-id",
            real_sample_files_import_to_loculus["batch_id"],
            "--timeline-file",
            str(real_sample_files_import_to_loculus["timeline_file"]),
            "--primer-file",
            str(real_sample_files_import_to_loculus["primer_file"]),
            "--output-fp",
            str(real_sample_files_import_to_loculus["output_file"]),
            "--reference",
            real_sample_files_import_to_loculus["reference"],
            "--skip-merge",
        ],
    )
    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout
    # The "Skip read pair merging: True" message is written to the logs, not stdout,
    # so we just verify the command was successful
