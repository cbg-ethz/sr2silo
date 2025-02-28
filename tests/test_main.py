from __future__ import annotations

from pathlib import Path

import pytest
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


@pytest.fixture
def sample_files(tmp_path):
    """Create sample files needed for testing."""
    input_file = tmp_path / "input.bam"
    input_file.write_text("")

    timeline_file = tmp_path / "timeline.tsv"
    timeline_file.write_text("sample_id\tdate\n1\t2024-01-01")

    primer_file = tmp_path / "primers.yaml"
    primer_file.write_text("primers:\n  - name: test\n    sequence: ACGT")

    output_file = tmp_path / "output.ndjson"

    return {
        "input_file": input_file,
        "timeline_file": timeline_file,
        "primer_file": primer_file,
        "output_file": output_file,
    }


def test_import_to_loculus_with_args(sample_files):
    """Test import-to-loculus with all required arguments."""
    result = runner.invoke(
        app,
        [
            "import-to-loculus",
            "--input-file",
            str(sample_files["input_file"]),
            "--sample-id",
            "test_sample",
            "--batch-id",
            "test_batch",
            "--timeline-file",
            str(sample_files["timeline_file"]),
            "--primer-file",
            str(sample_files["primer_file"]),
            "--output-fp",
            str(sample_files["output_file"]),
            "--reference",
            "sars-cov-2",
            "--no-upload",
        ],
    )

    # Command should start executing but might fail due to missing actual
    #  BAM file content
    # We're just testing the CLI interface here, not the actual processing
    assert "Starting V-PIPE to SILO conversion" in result.stdout


def test_import_to_loculus_with_upload_flag(sample_files):
    """Test import-to-loculus with upload flag."""
    result = runner.invoke(
        app,
        [
            "import-to-loculus",
            "--input-file",
            str(sample_files["input_file"]),
            "--sample-id",
            "test_sample",
            "--batch-id",
            "test_batch",
            "--timeline-file",
            str(sample_files["timeline_file"]),
            "--primer-file",
            str(sample_files["primer_file"]),
            "--output-fp",
            str(sample_files["output_file"]),
            "--upload",
        ],
    )

    assert "Starting V-PIPE to SILO conversion" in result.stdout
    # The command will proceed but actual upload will fail due to mock files


@pytest.fixture
def real_sample_files(tmp_path):
    """Get real sample files from the test data directory."""
    return {
        "input_file": Path(
            "./tests/data/samples/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam"
        ),
        "timeline_file": Path("./tests/data/samples/timeline_A1_05_2024_10_08.tsv"),
        "primer_file": Path("./tests/data/samples_large/primers.yaml"),
        "output_file": tmp_path / "silo_input.ndjson.zst",
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",
        "reference": "sars-cov-2",
    }


def test_import_to_loculus_with_real_files(real_sample_files):
    """Test import-to-loculus with real sample files."""

    result = runner.invoke(
        app,
        [
            "import-to-loculus",
            "--input-file",
            str(real_sample_files["input_file"]),
            "--sample-id",
            real_sample_files["sample_id"],
            "--batch-id",
            real_sample_files["batch_id"],
            "--timeline-file",
            str(real_sample_files["timeline_file"]),
            "--primer-file",
            str(real_sample_files["primer_file"]),
            "--output-fp",
            str(real_sample_files["output_file"]),
            "--reference",
            real_sample_files["reference"],
            "--no-upload",
        ],
    )

    print(result.stdout)
    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout


def test_import_to_loculus_with_real_files_and_upload(real_sample_files):
    """Test import-to-loculus with real sample files and upload flag."""
    result = runner.invoke(
        app,
        [
            "import-to-loculus",
            "--input-file",
            str(real_sample_files["input_file"]),
            "--sample-id",
            real_sample_files["sample_id"],
            "--batch-id",
            real_sample_files["batch_id"],
            "--timeline-file",
            str(real_sample_files["timeline_file"]),
            "--primer-file",
            str(real_sample_files["primer_file"]),
            "--output-fp",
            str(real_sample_files["output_file"]),
            "--reference",
            real_sample_files["reference"],
            "--upload",
        ],
    )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout
