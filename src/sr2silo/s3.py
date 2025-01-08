"""Functionality to interact with AWS S3."""

from __future__ import annotations

import bz2
import logging
import os
import shutil
from pathlib import Path

import boto3
from botocore.exceptions import NoCredentialsError
from moto import mock_aws

# Check if running in CI environment
is_CI = os.getenv("CI")


def compress_bz2(input_fp: Path, output_fp: Path) -> None:
    """Compress a file using BZ2 compression.

    Args:
        input_fp (Path): Path to the input file.
        output_fp (Path): Path to the output compressed file.
    """
    with open(input_fp, "rb") as f_in:
        with bz2.BZ2File(output_fp, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def get_aws_credentials():
    """Get AWS credentials from Docker secrets.

    Returns:
        Tuple[str, str, str]: AWS access key ID, AWS secret access key,
        and AWS default region.

    Raises:
        RuntimeError: If any of the required secrets are missing.
    """
    try:
        with open("/run/secrets/aws_access_key_id") as f:
            aws_access_key_id = f.read().strip()
        with open("/run/secrets/aws_secret_access_key") as f:
            aws_secret_access_key = f.read().strip()
        with open("/run/secrets/aws_default_region") as f:
            aws_default_region = f.read().strip()
    except FileNotFoundError as e:
        raise RuntimeError("Required secret is missing: " + str(e))

    return aws_access_key_id, aws_secret_access_key, aws_default_region


def get_s3_client():
    """Get an S3 client using AWS credentials from environment variables."""

    # Get AWS credentials from environment variables
    aws_access_key_id, aws_secret_access_key, aws_default_region = get_aws_credentials()

    # Create an S3 client
    s3_client = boto3.client(
        "s3",
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_default_region,
    )

    return s3_client


def upload_file_to_s3(file_name, bucket, object_name=None, client=None) -> bool:
    """Upload a file to an S3 bucket

    Args:
        file_name (str): Path to the file to upload.
        bucket (str): Bucket to upload to.
        object_name (str, optional): S3 object name.

    Returns:
        bool: True if the file was uploaded successfully, False otherwise.
    """
    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = file_name

    # If running in CI, mock the S3 upload
    if is_CI:
        logging.info("Running in CI environment, mocking S3 upload with moto.")
        with mock_aws():
            s3_client = boto3.client("s3", region_name="us-east-1")
            s3_client.create_bucket(Bucket=bucket)
            s3_client.upload_file(file_name, bucket, object_name or file_name)
            return True

    # If client was given, use it; otherwise, get the s3 client
    s3_client = client if client else get_s3_client()

    try:
        s3_client.upload_file(file_name, bucket, object_name)
    except NoCredentialsError:
        print("Credentials not available")
        return False
    return True


def download_file_from_s3(bucket_name, object_name, file_name=None):
    """Download a file from an S3 bucket.

    Args:
        bucket_name (str): Bucket to download from.
        object_name (str): S3 object name.
        file_name (str, optional): Path to the file to download to.
                            If not specified then object_name is used.
    """
    if file_name is None:
        file_name = object_name

    # get the s3 client
    s3_client = get_s3_client()

    try:
        s3_client.download_file(bucket_name, object_name, file_name)
        print(f"File {object_name} downloaded from {bucket_name} to {file_name}")
    except Exception as e:
        print(f"Error downloading file from S3: {e}")


def list_files_in_bucket(bucket_name):
    """List all files in an S3 bucket.

    Args:
        bucket_name (str): The name of the S3 bucket.
    """

    # get the s3 client
    s3_client = get_s3_client()

    try:
        response = s3_client.list_objects_v2(Bucket=bucket_name)
        if "Contents" in response:
            files = [obj["Key"] for obj in response["Contents"]]
            return files
        else:
            return []
    except Exception as e:
        print(f"Error listing files in bucket {bucket_name}: {e}")
        return []
