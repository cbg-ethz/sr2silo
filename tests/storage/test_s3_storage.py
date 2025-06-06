"""Testing the S3 module."""

from __future__ import annotations

import boto3
import pytest
from moto import mock_aws

from sr2silo.storage.s3_storage import upload_file_to_s3


@pytest.fixture
def s3_client():
    """Mocked S3 client."""
    with mock_aws():
        # Initialize the mock S3 service
        s3 = boto3.client("s3", region_name="us-east-1")
        # Create a mock bucket
        s3.create_bucket(Bucket="test-bucket")
        yield s3


def test_upload_file_to_s3(tmp_path, s3_client):
    """Test the upload_file_to_s3 function."""
    # Create a test file
    test_file = tmp_path / "test.txt"
    test_file.write_text("This is a test file.")

    # Upload the test file to the mock S3 bucket
    upload_file_to_s3(str(test_file), "test-bucket", "test.txt", client=s3_client)

    # Check that the file exists in the mock S3 bucket
    response = s3_client.list_objects_v2(Bucket="test-bucket")
    assert "Contents" in response
    assert any(obj["Key"] == "test.txt" for obj in response["Contents"])
