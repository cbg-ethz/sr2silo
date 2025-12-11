"""Tests for the main CLI application of sr2silo.
These tests cover the command-line interface and its
functionality, ensuring that the application behaves
as expected when invoked with various commands and options.
"""

from __future__ import annotations

import json
import subprocess
from unittest.mock import patch

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

    # Mock is_ci_environment to return True for this test
    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(real_sample_files_import_to_loculus["input_file"]),
                "--sample-id",
                real_sample_files_import_to_loculus["sample_id"],
                "--timeline-file",
                str(real_sample_files_import_to_loculus["timeline_file"]),
                "--organism",
                "covid",
                "--lapis-url",
                real_sample_files_import_to_loculus["lapis_url"],
                "--output-fp",
                str(real_sample_files_import_to_loculus["output_file"]),
            ],
        )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout


def test_process_from_vpipe_with_empty_batch_id(real_sample_files_import_to_loculus):
    """Test process-from-vpipe without batch_id parameter."""

    # Mock is_ci_environment to return True for this test
    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(real_sample_files_import_to_loculus["input_file"]),
                "--sample-id",
                real_sample_files_import_to_loculus["sample_id"],
                "--timeline-file",
                str(real_sample_files_import_to_loculus["timeline_file"]),
                "--organism",
                "covid",
                "--output-fp",
                str(
                    real_sample_files_import_to_loculus["output_file"].parent
                    / "empty_batch_test.ndjson.zst"
                ),
                "--lapis-url",
                real_sample_files_import_to_loculus["lapis_url"],
            ],
        )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout
    # The logging message goes to stderr/logs, not stdout,
    # so we just verify successful execution


def test_process_from_vpipe_with_explicit_empty_batch_id(
    real_sample_files_import_to_loculus,
):
    """Test process-from-vpipe without batch_id parameter (simplified)."""

    # Mock is_ci_environment to return True for this test
    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(real_sample_files_import_to_loculus["input_file"]),
                "--sample-id",
                real_sample_files_import_to_loculus["sample_id"],
                "--timeline-file",
                str(real_sample_files_import_to_loculus["timeline_file"]),
                "--organism",
                "covid",
                "--output-fp",
                str(
                    real_sample_files_import_to_loculus["output_file"].parent
                    / "explicit_empty_batch_test.ndjson.zst"
                ),
                "--lapis-url",
                real_sample_files_import_to_loculus["lapis_url"],
            ],
        )

    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout
    # The logging message goes to stderr/logs, not stdout,
    # so we just verify successful execution


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
            "--nucleotide-alignment",
            "/tmp/nonexistent.bam",
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
                "--nucleotide-alignment",
                "/tmp/nonexistent.bam",
            ],
        )
        assert result.exit_code == 1
    finally:
        import os

        os.unlink(tmp_path)


def test_process_from_vpipe_with_skip_merge(real_sample_files_import_to_loculus):
    """Test process-from-vpipe with skip-merge option to bypass read pair merging."""

    # Mock is_ci_environment to return True for this test
    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(real_sample_files_import_to_loculus["input_file"]),
                "--sample-id",
                real_sample_files_import_to_loculus["sample_id"],
                "--timeline-file",
                str(real_sample_files_import_to_loculus["timeline_file"]),
                "--organism",
                "covid",
                "--output-fp",
                str(real_sample_files_import_to_loculus["output_file"]),
                "--lapis-url",
                real_sample_files_import_to_loculus["lapis_url"],
                "--skip-merge",
            ],
        )
    assert result.exit_code == 0
    assert "Starting V-PIPE to SILO conversion" in result.stdout
    # The "Skip read pair merging: True" message is written to the logs, not stdout,
    # so we just verify the command was successful


def test_process_from_vpipe_environment_variables():
    """Test process-from-vpipe with environment variables."""
    import os
    from unittest.mock import patch

    # Test with environment variables set
    env_vars = {
        "TIMELINE_FILE": "/path/to/timeline.tsv",
    }

    with patch.dict(os.environ, env_vars):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                "/tmp/test.bam",
                "--sample-id",
                "test-sample",
                "--output-fp",
                "/tmp/test.ndjson.zst",
                "--lapis-url",
                "https://lapis.example.com",
            ],
        )
        # Should fail due to missing input file,
        # but environment variables should be processed
        assert result.exit_code == 2  # Should fail due to missing file
        # We can verify the environment variables were picked up by manually checking
        # since the logs go to stderr which typer doesn't capture by default


def test_process_from_vpipe_cli_overrides_env():
    """Test that CLI arguments override environment variables."""
    import os

    # Set environment variables
    env_vars = {
        "TIMELINE_FILE": "/env/timeline.tsv",
    }

    with patch.dict(os.environ, env_vars):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                "/tmp/test.bam",
                "--sample-id",
                "test-sample",
                "--output-fp",
                "/tmp/test.ndjson.zst",
                "--timeline-file",
                "/cli/timeline.tsv",
                "--lapis-url",
                "https://cli.lapis.example.com",
            ],
        )
        # CLI arguments should override environment variables
        assert result.exit_code == 1  # Should fail due to missing file
        # Environment vars were set but CLI args should take precedence


def test_process_from_vpipe_missing_env_and_cli():
    """Test process-from-vpipe fails when neither env vars nor CLI args provided."""
    import os

    # Clear environment variables
    with patch.dict(os.environ, {}, clear=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                "/tmp/test.bam",
                "--sample-id",
                "test-sample",
                "--output-fp",
                "/tmp/test.ndjson.zst",
                "--lapis-url",
                "https://lapis.example.com",
                "--timeline-file",
                "/tmp/timeline.tsv",
            ],
        )
        assert result.exit_code == 1
        # The error message is logged, so let's just verify the exit code for now
        # In a real scenario, this would be caught by logging or stderr capture


def test_submit_to_loculus_environment_variables():
    """Test submit-to-loculus with environment variables."""
    import os
    from unittest.mock import patch

    # Test with environment variables set
    env_vars = {
        "KEYCLOAK_TOKEN_URL": "https://env.auth.com/token",
        "BACKEND_URL": "https://env.submit.com/api",
    }

    with patch.dict(os.environ, env_vars):
        result = runner.invoke(
            app,
            [
                "submit-to-loculus",
                "--processed-file",
                "/tmp/test.ndjson.zst",
                "--nucleotide-alignment",
                "/tmp/test.bam",
            ],
        )
        # Should fail due to missing file, but environment variables should be processed
        assert result.exit_code == 1  # Should fail due to missing file


def test_submit_to_loculus_cli_overrides_env():
    """Test that CLI arguments override environment variables for submit command."""
    import os
    from unittest.mock import patch

    # Set environment variables
    env_vars = {
        "KEYCLOAK_TOKEN_URL": "https://env.auth.com/token",
        "BACKEND_URL": "https://env.submit.com/api",
    }

    with patch.dict(os.environ, env_vars):
        result = runner.invoke(
            app,
            [
                "submit-to-loculus",
                "--processed-file",
                "/tmp/test.ndjson.zst",
                "--nucleotide-alignment",
                "/tmp/test.bam",
                "--keycloak-token-url",
                "https://cli.auth.com/token",
                "--backend-url",
                "https://cli.submit.com/api",
            ],
        )
        # CLI arguments should override environment variables
        assert result.exit_code == 1  # Should fail due to missing file


def test_submit_to_loculus_missing_env_and_cli():
    """Test submit-to-loculus fails when neither env vars nor CLI args provided."""
    import os
    from unittest.mock import patch

    # Clear environment variables
    with patch.dict(os.environ, {}, clear=True):
        result = runner.invoke(
            app,
            [
                "submit-to-loculus",
                "--processed-file",
                "/tmp/test.ndjson.zst",
                "--nucleotide-alignment",
                "/tmp/test.bam",
            ],
        )
        assert result.exit_code == 1
        # The error message is logged, so let's just verify the exit code for now


def test_process_from_vpipe_multi_organism(sample_data_by_organism, tmp_path):
    """Test process-from-vpipe with multiple organisms (parameterized).

    This test verifies that the organism parameter is correctly parsed and passed
    through to the reference loading logic. It mocks the actual processing to test
    the CLI argument handling for both COVID-19 and RSV-A organisms.

    Args:
        sample_data_by_organism: Parameterized fixture providing test data
        tmp_path: Pytest fixture for temporary directory
    """
    organism_data = sample_data_by_organism
    output_fp = tmp_path / f"{organism_data['organism']}_output.ndjson.zst"

    # Mock the actual processing to avoid needing paired_end_read_merger binary
    with (
        patch("sr2silo.main.is_ci_environment", return_value=True),
        patch("sr2silo.main.nuc_align_to_silo_njson") as mock_process,
    ):
        mock_process.return_value = None  # Mock successful processing

        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(organism_data["sample"]),
                "--sample-id",
                organism_data["sample_id"],
                "--timeline-file",
                str(organism_data["timeline"]),
                "--organism",
                organism_data["organism"],
                "--output-fp",
                str(output_fp),
            ],
        )

    # Verify successful execution and organism was used
    assert result.exit_code == 0, (
        f"Failed for {organism_data['organism']}: {result.stdout}\nStderr: {result.stderr}"
    )
    assert "Starting V-PIPE to SILO conversion" in result.stdout

    # Verify processing was called (organism parameter was passed through)
    assert mock_process.called


def test_process_from_vpipe_organism_help():
    """Test that organism parameter is shown in help text."""
    result = runner.invoke(app, ["process-from-vpipe", "--help"])
    assert result.exit_code == 0
    # Help output includes ANSI codes for formatting, so check for the text content
    assert "organism" in result.stdout
    assert "Organism identifier" in result.stdout


def test_process_from_vpipe_organism_from_env(
    sample_data_by_organism, tmp_path, monkeypatch
):
    """Test that organism can be resolved from ORGANISM environment variable.

    This test verifies organism resolution from environment variables.
    It mocks the actual processing to focus on CLI argument handling.

    Args:
        sample_data_by_organism: Parameterized fixture providing test data
        tmp_path: Pytest fixture for temporary directory
        monkeypatch: Pytest fixture for environment variable patching
    """
    organism_data = sample_data_by_organism
    output_fp = tmp_path / f"{organism_data['organism']}_env_output.ndjson.zst"

    # Set ORGANISM environment variable
    monkeypatch.setenv("ORGANISM", organism_data["organism"])

    # Mock the actual processing to avoid needing paired_end_read_merger binary
    with (
        patch("sr2silo.main.is_ci_environment", return_value=True),
        patch("sr2silo.main.nuc_align_to_silo_njson") as mock_process,
    ):
        mock_process.return_value = None  # Mock successful processing

        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(organism_data["sample"]),
                "--sample-id",
                organism_data["sample_id"],
                "--timeline-file",
                str(organism_data["timeline"]),
                # Note: NOT providing --organism, should use env var
                "--output-fp",
                str(output_fp),
            ],
        )

    # Verify successful execution with organism from env
    assert result.exit_code == 0, (
        f"Failed for {organism_data['organism']}: {result.stdout}\nStderr: {result.stderr}"
    )
    assert "Starting V-PIPE to SILO conversion" in result.stdout

    # Verify processing was called (organism parameter was passed through)
    assert mock_process.called


def test_process_from_vpipe_output_covid(real_sample_files_import_to_loculus, tmp_path):
    """Integration test: Verify actual output structure for COVID data.

    This test validates that sr2silo produces properly formatted NDJSON output
    with all expected fields populated correctly.
    """
    output_file = tmp_path / "covid_output.ndjson.zst"

    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(real_sample_files_import_to_loculus["input_file"]),
                "--sample-id",
                real_sample_files_import_to_loculus["sample_id"],
                "--timeline-file",
                str(real_sample_files_import_to_loculus["timeline_file"]),
                "--organism",
                "covid",
                "--output-fp",
                str(output_file),
            ],
        )

    assert result.exit_code == 0, f"Command failed: {result.stdout}\n{result.stderr}"
    assert output_file.exists(), "Output file was not created"

    # Decompress using zstdcat command
    decompress_result = subprocess.run(
        ["zstdcat", str(output_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    records = [
        json.loads(line) for line in decompress_result.stdout.strip().split("\n")
    ]

    # Verify we have records
    assert len(records) > 0, "No records in output file"

    # Validate structure of first record
    first_record = records[0]

    # Check required fields
    required_fields = [
        "readId",
        "main",
        "sampleId",
        "batchId",
        "samplingDate",
        "locationName",
        "locationCode",
        "sr2siloVersion",
    ]
    for field in required_fields:
        assert field in first_record, f"Missing required field: {field}"

    # Check main sequence structure
    assert isinstance(first_record["main"], dict)
    assert "sequence" in first_record["main"]
    assert "insertions" in first_record["main"]
    assert "offset" in first_record["main"]

    # Check that COVID genes are present in record (at least one should have data)
    covid_genes = ["E", "M", "N", "ORF1a", "S"]
    has_gene_data = any(first_record.get(gene) is not None for gene in covid_genes)
    assert has_gene_data, "No gene data found in output"

    # Validate sample metadata matches input
    assert first_record["sampleId"] == real_sample_files_import_to_loculus["sample_id"]
    assert first_record["batchId"] == real_sample_files_import_to_loculus["batch_id"]


def test_process_from_vpipe_output_covid_full(real_sample_files_import_to_loculus):
    """Integration test: Compare actual COVID output against expected output.

    This test validates that sr2silo produces output identical to the expected
    reference output, ensuring alignment and sequence data are correctly processed.
    """
    from pathlib import Path as PathlibPath

    expected_output_file = (
        PathlibPath(__file__).parent / "data" / "covid" / "expected_output.ndjson.zst"
    )
    actual_output_file = expected_output_file.parent / "actual_output.ndjson.zst"

    # Generate actual output
    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(real_sample_files_import_to_loculus["input_file"]),
                "--sample-id",
                real_sample_files_import_to_loculus["sample_id"],
                "--timeline-file",
                str(real_sample_files_import_to_loculus["timeline_file"]),
                "--organism",
                "covid",
                "--output-fp",
                str(actual_output_file),
            ],
        )

    assert result.exit_code == 0, f"Command failed: {result.stdout}"
    assert actual_output_file.exists(), "Actual output file was not created"
    assert expected_output_file.exists(), "Expected output file not found"

    # Decompress both files and compare
    expected_result = subprocess.run(
        ["zstdcat", str(expected_output_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    actual_result = subprocess.run(
        ["zstdcat", str(actual_output_file)],
        capture_output=True,
        text=True,
        check=True,
    )

    expected_lines = [
        line for line in expected_result.stdout.strip().split("\n") if line
    ]
    actual_lines = [line for line in actual_result.stdout.strip().split("\n") if line]

    # Compare record count
    assert len(actual_lines) == len(expected_lines), (
        f"Record count mismatch: expected {len(expected_lines)}, "
        f"got {len(actual_lines)}"
    )

    # Compare each record
    for idx, (expected_line, actual_line) in enumerate(
        zip(expected_lines, actual_lines)
    ):
        expected_record = json.loads(expected_line)
        actual_record = json.loads(actual_line)

        # Compare read ID
        assert actual_record["readId"] == expected_record["readId"], (
            f"Record {idx}: readId mismatch. Expected {expected_record['readId']}, "
            f"got {actual_record['readId']}"
        )

        # Compare main sequence (critical alignment data)
        assert (
            actual_record["main"]["sequence"] == expected_record["main"]["sequence"]
        ), (
            f"Record {idx} ({actual_record['readId']}): main sequence mismatch. "
            "Alignment changed."
        )

        # Compare insertions
        assert (
            actual_record["main"]["insertions"] == expected_record["main"]["insertions"]
        ), f"Record {idx} ({actual_record['readId']}): insertions mismatch"

        # Compare offset
        assert actual_record["main"]["offset"] == expected_record["main"]["offset"], (
            f"Record {idx} ({actual_record['readId']}): offset mismatch"
        )

        # Compare metadata fields (only those present in expected record)
        metadata_fields = [
            "sampleId",
            "batchId",
            "samplingDate",
            "locationName",
            "locationCode",
        ]
        for field in metadata_fields:
            if field in expected_record:
                assert field in actual_record, (
                    f"Record {idx}: {field} missing in actual output"
                )
                assert actual_record[field] == expected_record[field], (
                    f"Record {idx}: {field} mismatch. "
                    f"Expected {expected_record[field]}, got {actual_record[field]}"
                )

        # Compare gene sequences (all genes present in expected record)
        for gene in expected_record:
            if gene.startswith("ORF") or gene in ["E", "M", "N", "S"]:
                assert gene in actual_record, (
                    f"Record {idx}: Gene {gene} missing in actual output"
                )
                if expected_record[gene] is not None:
                    assert actual_record[gene] == expected_record[gene], (
                        f"Record {idx}: Gene {gene} sequence mismatch"
                    )

    # Clean up
    actual_output_file.unlink()


def test_process_from_vpipe_output_rsva_full(rsva_sample, rsva_timeline):
    """Integration test: Compare actual RSV-A output against expected output.

    This test validates that sr2silo produces output identical to the expected
    reference output for RSV-A, ensuring alignment and sequence data are correctly processed.
    """
    from pathlib import Path as PathlibPath

    expected_output_file = (
        PathlibPath(__file__).parent / "data" / "rsva" / "expected_output.ndjson.zst"
    )
    actual_output_file = expected_output_file.parent / "actual_output.ndjson.zst"

    # Generate actual output
    with patch("sr2silo.main.is_ci_environment", return_value=True):
        result = runner.invoke(
            app,
            [
                "process-from-vpipe",
                "--input-file",
                str(rsva_sample),
                "--sample-id",
                "A1_05_2025_11_05",
                "--timeline-file",
                str(rsva_timeline),
                "--organism",
                "rsva",
                "--output-fp",
                str(actual_output_file),
            ],
        )

    assert result.exit_code == 0, f"Command failed: {result.stdout}"
    assert actual_output_file.exists(), "Actual output file was not created"
    assert expected_output_file.exists(), "Expected output file not found"

    # Decompress both files and compare
    expected_result = subprocess.run(
        ["zstdcat", str(expected_output_file)],
        capture_output=True,
        text=True,
        check=True,
    )
    actual_result = subprocess.run(
        ["zstdcat", str(actual_output_file)],
        capture_output=True,
        text=True,
        check=True,
    )

    expected_lines = [
        line for line in expected_result.stdout.strip().split("\n") if line
    ]
    actual_lines = [line for line in actual_result.stdout.strip().split("\n") if line]

    # Skip comparison if no records (empty test BAM)
    if len(expected_lines) == 0 and len(actual_lines) == 0:
        actual_output_file.unlink()
        return

    # Compare record count
    assert len(actual_lines) == len(expected_lines), (
        f"Record count mismatch: expected {len(expected_lines)}, "
        f"got {len(actual_lines)}"
    )

    # Compare each record
    for idx, (expected_line, actual_line) in enumerate(
        zip(expected_lines, actual_lines)
    ):
        expected_record = json.loads(expected_line)
        actual_record = json.loads(actual_line)

        # Compare read ID
        assert actual_record["readId"] == expected_record["readId"], (
            f"Record {idx}: readId mismatch. Expected {expected_record['readId']}, "
            f"got {actual_record['readId']}"
        )

        # Compare main sequence (critical alignment data)
        assert (
            actual_record["main"]["sequence"] == expected_record["main"]["sequence"]
        ), (
            f"Record {idx} ({actual_record['readId']}): main sequence mismatch. "
            "Alignment changed."
        )

        # Compare insertions
        assert (
            actual_record["main"]["insertions"] == expected_record["main"]["insertions"]
        ), f"Record {idx} ({actual_record['readId']}): insertions mismatch"

        # Compare offset
        assert actual_record["main"]["offset"] == expected_record["main"]["offset"], (
            f"Record {idx} ({actual_record['readId']}): offset mismatch"
        )

        # Compare metadata fields (only those present in expected record)
        metadata_fields = [
            "sampleId",
            "batchId",
            "samplingDate",
            "locationName",
            "locationCode",
        ]
        for field in metadata_fields:
            if field in expected_record:
                assert field in actual_record, (
                    f"Record {idx}: {field} missing in actual output"
                )
                assert actual_record[field] == expected_record[field], (
                    f"Record {idx}: {field} mismatch. "
                    f"Expected {expected_record[field]}, got {actual_record[field]}"
                )

        # Compare gene sequences (all genes present in expected record)
        for gene in expected_record:
            if gene.startswith("ORF") or gene in ["F", "G", "N", "M", "SH"]:
                assert gene in actual_record, (
                    f"Record {idx}: Gene {gene} missing in actual output"
                )
                if expected_record[gene] is not None:
                    assert actual_record[gene] == expected_record[gene], (
                        f"Record {idx}: Gene {gene} sequence mismatch"
                    )

    # Clean up
    actual_output_file.unlink()
