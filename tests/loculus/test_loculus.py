"""Tests for loculus module functions."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml

from sr2silo.loculus.loculus import LoculusClient, Submission


@pytest.fixture
def test_silo_input_uncompressed():
    """Fixture providing path to uncompressed test SILO input file."""
    return Path(__file__).parent / "test_silo_input.ndjson"


@pytest.fixture
def test_silo_input_compressed():
    """Fixture providing path to compressed test SILO input file."""
    return Path(__file__).parent / "test_silo_input.ndjson.zst"


def test_parse_metadata_uncompressed(test_silo_input_uncompressed):
    """Test parsing metadata from uncompressed SILO input file."""
    metadata = Submission.parse_metadata(test_silo_input_uncompressed)

    # Check that metadata is returned as dictionary
    assert isinstance(metadata, dict)

    # Check expected metadata fields in snake_case (excluding read_id)
    expected_fields = {
        "sample_id": "A1_05_2025_06_18",
        "batch_id": "20250711_2443602573",
        "location_code": "5",
        "location_name": "Lugano (TI)",
        "sampling_date": "2025-06-18",
        "sr2silo_version": "1.3.0 (v1.1.0-5-g31b0623)",
    }

    for key, expected_value in expected_fields.items():
        assert key in metadata
        assert metadata[key] == expected_value

    # Ensure primary key field is excluded
    assert "read_id" not in metadata


def test_parse_metadata_compressed(test_silo_input_compressed):
    """Test parsing metadata from compressed SILO input file."""
    metadata = Submission.parse_metadata(test_silo_input_compressed)

    # Check that metadata is returned as dictionary
    assert isinstance(metadata, dict)

    # Check expected metadata fields in snake_case (excluding read_id)
    expected_fields = {
        "sample_id": "A1_05_2025_06_18",
        "batch_id": "20250711_2443602573",
        "location_code": "5",
        "location_name": "Lugano (TI)",
        "sampling_date": "2025-06-18",
        "sr2silo_version": "1.3.0 (v1.1.0-5-g31b0623)",
    }

    for key, expected_value in expected_fields.items():
        assert key in metadata
        assert metadata[key] == expected_value

    # Ensure primary key field is excluded
    assert "read_id" not in metadata


def test_count_reads_uncompressed(test_silo_input_uncompressed):
    """Test counting reads in uncompressed SILO input file."""
    count = Submission.count_reads(test_silo_input_uncompressed)

    # The test file has 2 lines (no final newline), so expect 2 reads
    assert count == 3


def test_count_reads_compressed(test_silo_input_compressed):
    """Test counting reads in compressed SILO input file."""
    count = Submission.count_reads(test_silo_input_compressed)

    # The test file has 3 lines, so expect 3 reads
    assert count == 3


def test_create_metadata_file(test_silo_input_uncompressed):
    """Test creating metadata TSV file."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy test file to temp directory to simulate processed file
        temp_processed_file = Path(temp_dir) / "test_processed.ndjson"
        temp_processed_file.write_text(test_silo_input_uncompressed.read_text())

        # Test without count_reads
        metadata_file, submission_id = Submission.create_metadata_file(
            temp_processed_file
        )

        # Check that file was created
        assert metadata_file.exists()
        assert metadata_file.name.startswith("metadata_")
        assert metadata_file.suffix == ".tsv"

        # Check that submission directory was created
        assert metadata_file.parent.name == "submission"

        # Check file contents
        with open(metadata_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            fieldnames = reader.fieldnames or []

            # Should have exactly 1 row of data (wide format)
            assert len(rows) == 1
            assert "submissionId" in fieldnames
            assert "date" in fieldnames

            # Check submission ID matches
            first_row = rows[0]
            assert first_row["submissionId"] == submission_id
            assert first_row["date"]  # Should have a date

            # Check that metadata fields are present as columns
            expected_metadata_fields = {
                "sampleId",
                "batchId",
                "locationCode",
                "locationName",
                "samplingDate",
                "sr2siloVersion",
            }
            for field in expected_metadata_fields:
                assert field in fieldnames, (
                    f"Expected metadata field '{field}' not found in columns"
                )
                assert first_row[field]  # Should have a value

        # Test with count_reads=True
        metadata_file2, submission_id2 = Submission.create_metadata_file(
            temp_processed_file, count_reads=True
        )

        # Should create a different file with different submission ID
        assert metadata_file2 != metadata_file
        assert submission_id2 != submission_id
        assert metadata_file2.exists()

        # Check that countReads column is included when count_reads=True
        with open(metadata_file2, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            fieldnames = reader.fieldnames or []

            assert "countSiloReads" in fieldnames
            first_row = rows[0]
            assert first_row["countSiloReads"]  # Should have a count value
            assert int(first_row["countSiloReads"]) > 0  # Should be a positive number


def test_metadata_fields_match_organism_config(test_silo_input_uncompressed):
    """Test that created metadata file fields match organism_conf.yml schema."""
    # Load organism configuration
    organism_conf_path = (
        Path(__file__).parent.parent.parent
        / "resources"
        / "loculus"
        / "organism_conf.yml"
    )
    assert organism_conf_path.exists(), (
        f"organism_conf.yml not found at {organism_conf_path}"
    )

    with open(organism_conf_path) as f:
        org_config = yaml.safe_load(f)

    # Extract metadata fields from organism config
    covid_schema = org_config["organisms"]["covid"]["schema"]
    config_metadata_fields = {field["name"] for field in covid_schema["metadata"]}

    # Remove primerProtocol from config fields as we no longer include it
    config_metadata_fields.discard("primerProtocol")

    # Create metadata file
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_processed_file = Path(temp_dir) / "test_processed.ndjson"
        temp_processed_file.write_text(test_silo_input_uncompressed.read_text())

        metadata_file, _ = Submission.create_metadata_file(
            temp_processed_file, count_reads=True
        )

        # Read the created metadata file
        with open(metadata_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = set(reader.fieldnames or [])

        # Extract metadata field names (excluding submissionId which is a special field)
        created_metadata_fields = fieldnames - {"submissionId"}

        # Check that all config fields are present in created file
        missing_fields = config_metadata_fields - created_metadata_fields
        assert not missing_fields, f"Missing fields in metadata file: {missing_fields}"

        # Check that no unexpected fields are present
        unexpected_fields = created_metadata_fields - config_metadata_fields
        assert not unexpected_fields, (
            f"Unexpected fields in metadata file: {unexpected_fields}"
        )

        # Verify exact match
        assert created_metadata_fields == config_metadata_fields, (
            f"Metadata fields don't match. Expected: {config_metadata_fields}, Got: {created_metadata_fields}"
        )


class TestLoculusClientApproveAll:
    """Tests for the LoculusClient.approve_all() method."""

    def test_approve_all_requires_authentication(self):
        """Test that approve_all raises exception when not authenticated."""
        client = LoculusClient(
            token_url="https://example.com/token",
            backend_url="https://example.com/backend",
            organism="covid",
        )
        # Token is None by default
        with pytest.raises(Exception) as exc_info:
            client.approve_all(username="testuser")
        assert "Authentication required" in str(exc_info.value)

    @patch("sr2silo.loculus.loculus.is_ci_environment")
    def test_approve_all_ci_environment(self, mock_ci_env):
        """Test that approve_all returns mock response in CI environment."""
        mock_ci_env.return_value = True

        client = LoculusClient(
            token_url="https://example.com/token",
            backend_url="https://example.com/backend",
            organism="covid",
        )
        client.token = "dummy_token"

        result = client.approve_all(username="testuser")

        assert result["status"] == "success"
        assert "Mock approval" in result["message"]

    @patch("sr2silo.loculus.loculus.requests.post")
    @patch("sr2silo.loculus.loculus.is_ci_environment")
    def test_approve_all_success(self, mock_ci_env, mock_post):
        """Test successful approval call."""
        mock_ci_env.return_value = False

        # Mock successful response
        mock_response = MagicMock()
        mock_response.json.return_value = {"approved": 5}
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response

        client = LoculusClient(
            token_url="https://example.com/token",
            backend_url="https://backend.example.com",
            organism="covid",
        )
        client.token = "test_token"

        result = client.approve_all(username="testuser")

        # Verify the API was called correctly
        mock_post.assert_called_once()
        call_args = mock_post.call_args

        # Check URL
        assert (
            call_args[0][0]
            == "https://backend.example.com/covid/approve-processed-data"
        )

        # Check headers
        headers = call_args[1]["headers"]
        assert headers["Authorization"] == "Bearer test_token"
        assert headers["Content-Type"] == "application/json"

        # Check payload
        payload = call_args[1]["json"]
        assert payload["scope"] == "ALL"
        assert payload["submitterNamesFilter"] == ["testuser"]

        # Check result
        assert result == {"approved": 5}

    @patch("sr2silo.loculus.loculus.requests.post")
    @patch("sr2silo.loculus.loculus.is_ci_environment")
    def test_approve_all_http_error(self, mock_ci_env, mock_post):
        """Test approve_all handles HTTP errors properly."""
        mock_ci_env.return_value = False

        # Mock HTTP error response
        mock_response = MagicMock()
        mock_response.status_code = 403
        mock_response.text = "Forbidden"
        mock_response.raise_for_status.side_effect = Exception("HTTP Error")
        mock_post.return_value = mock_response

        client = LoculusClient(
            token_url="https://example.com/token",
            backend_url="https://backend.example.com",
            organism="covid",
        )
        client.token = "test_token"

        with pytest.raises(Exception) as exc_info:
            client.approve_all(username="testuser")

        assert "Failed to approve" in str(exc_info.value) or "HTTP Error" in str(
            exc_info.value
        )


class TestSubmitWithAutoRelease:
    """Integration tests for submit() with auto_release enabled."""

    @patch("sr2silo.submit_to_loculus.time.sleep")
    @patch("sr2silo.loculus.loculus.is_ci_environment")
    @patch("sr2silo.submit_to_loculus.is_ci_environment")
    def test_submit_with_auto_release_calls_approve(
        self, mock_submit_ci, mock_loculus_ci, mock_sleep, test_silo_input_uncompressed
    ):
        """Test that submit() with auto_release=True calls approve_all after submission."""
        # Both CI environment checks need to return True
        mock_submit_ci.return_value = True
        mock_loculus_ci.return_value = True

        # Create a temporary processed file and nucleotide alignment
        with tempfile.TemporaryDirectory() as temp_dir:
            # Copy test file to temp directory
            temp_processed = Path(temp_dir) / "test.ndjson.zst"

            # Create a minimal zstd-compressed file
            import zstandard as zstd

            content = test_silo_input_uncompressed.read_text()
            cctx = zstd.ZstdCompressor()
            compressed = cctx.compress(content.encode())
            temp_processed.write_bytes(compressed)

            # Create a dummy BAM file
            temp_bam = Path(temp_dir) / "test.bam"
            temp_bam.write_bytes(b"dummy bam content")

            # Import and call submit
            from sr2silo.submit_to_loculus import submit

            result = submit(
                processed_file=temp_processed,
                nucleotide_alignment=temp_bam,
                keycloak_token_url="https://example.com/token",
                backend_url="https://example.com/backend",
                group_id=1,
                organism="covid",
                username="testuser",
                password="testpass",
                auto_release=True,
                release_delay=0,  # No delay for testing
            )

            # Verify submission succeeded
            assert result is True

            # Verify sleep was called (even with 0 delay)
            mock_sleep.assert_called_once_with(0)

    @patch("sr2silo.submit_to_loculus.time.sleep")
    @patch("sr2silo.loculus.loculus.is_ci_environment")
    @patch("sr2silo.submit_to_loculus.is_ci_environment")
    def test_submit_without_auto_release_skips_approve(
        self, mock_submit_ci, mock_loculus_ci, mock_sleep, test_silo_input_uncompressed
    ):
        """Test that submit() without auto_release does not call approve_all."""
        mock_submit_ci.return_value = True
        mock_loculus_ci.return_value = True

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_processed = Path(temp_dir) / "test.ndjson.zst"

            import zstandard as zstd

            content = test_silo_input_uncompressed.read_text()
            cctx = zstd.ZstdCompressor()
            compressed = cctx.compress(content.encode())
            temp_processed.write_bytes(compressed)

            temp_bam = Path(temp_dir) / "test.bam"
            temp_bam.write_bytes(b"dummy bam content")

            from sr2silo.submit_to_loculus import submit

            result = submit(
                processed_file=temp_processed,
                nucleotide_alignment=temp_bam,
                keycloak_token_url="https://example.com/token",
                backend_url="https://example.com/backend",
                group_id=1,
                organism="covid",
                username="testuser",
                password="testpass",
                auto_release=False,  # Disabled
                release_delay=180,
            )

            assert result is True

            # Verify sleep was NOT called (auto_release=False)
            mock_sleep.assert_not_called()
