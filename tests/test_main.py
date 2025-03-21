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


def test_import_to_loculus_help():
    """Test the help output for import-to-loculus command."""
    result = runner.invoke(app, ["import-to-loculus", "--help"])
    assert result.exit_code == 0
    assert "V-PIPE to SILO conversion" in result.stdout


def test_import_to_loculus_missing_required():
    """Test import-to-loculus fails when required arguments are missing."""
    result = runner.invoke(app, ["import-to-loculus"])
    assert result.exit_code != 0  # Should fail due to missing required options


def test_import_to_loculus_with_real_files(real_sample_files_import_to_loculus):
    """Test import-to-loculus with real sample files."""

    result = runner.invoke(
        app,
        [
            "import-to-loculus",
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
            "--no-upload",
        ],
    )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout


def test_import_to_loculus_with_real_files_and_upload(
    real_sample_files_import_to_loculus,
):
    """Test import-to-loculus with real sample files and upload flag in a CI
    environment.
    """
    env = os.environ.copy()
    env["CI"] = "true"

    result = runner.invoke(
        app,
        [
            "import-to-loculus",
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
            "--upload",
        ],
        env=env,
    )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout
