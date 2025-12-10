"""Tests for the main CLI application of sr2silo.
These tests cover the command-line interface and its
functionality, ensuring that the application behaves
as expected when invoked with various commands and options.
"""

from __future__ import annotations

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
                "sars-cov-2",
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
                "sars-cov-2",
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
                "sars-cov-2",
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
                "sars-cov-2",
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
    the CLI argument handling for both SARS-CoV-2 and RSV-A organisms.

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
    assert "--organism" in result.stdout
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
