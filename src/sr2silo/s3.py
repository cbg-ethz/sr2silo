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

from sr2silo.config import is_ci_environment

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

# Set logging level to WARNING for boto3 and botocore
logging.getLogger("boto3").setLevel(logging.WARNING)
logging.getLogger("botocore").setLevel(logging.WARNING)
logging.getLogger("s3transfer").setLevel(logging.WARNING)


def compress_bz2(input_fp: Path, output_fp: Path) -> None:
    """Compress a file using BZ2 compression.

    Args:
        input_fp (Path): Path to the input file.
        output_fp (Path): Path to the output compressed file.
    """
    with open(input_fp, "rb") as f_in:
        with bz2.BZ2File(output_fp, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def _get_aws_credentials(secrets_dir: Path) -> tuple[str, str, str]:
    """Get AWS credentials secrets.

    Returns:
        Tuple[str, str, str]: AWS access key ID, AWS secret access key,
        and AWS default region.

    Raises:
        RuntimeError: If any of the required secrets are missing.
    """
    try:
        with open(secrets_dir / "aws_access_key_id") as f:
            aws_access_key_id = f.read().strip()
        with open(secrets_dir / "aws_secret_access_key") as f:
            aws_secret_access_key = f.read().strip()
        with open(secrets_dir / "aws_default_region") as f:
            aws_default_region = f.read().strip()
    except FileNotFoundError as e:
        raise RuntimeError("Required secret is missing: " + str(e))

    return aws_access_key_id, aws_secret_access_key, aws_default_region


def get_aws_credentials_from_docker_secrets():
    """Get AWS credentials from Docker secrets.

    Returns:
        Tuple[str, str, str]: AWS access key ID, AWS secret access key,
        and AWS default region.

    Raises:
        RuntimeError: If any of the required secrets are missing.
    """
    return _get_aws_credentials(Path("/run/secrets/"))


def _is_running_in_docker():
    # Check if running in Docker by verifying the existence of /.dockerenv
    if os.path.exists("/.dockerenv"):
        return True
    # Alternatively, check /proc/1/cgroup for docker indication
    try:
        with open("/proc/1/cgroup", "rt") as f:
            if any("docker" in line for line in f):
                return True
    except Exception:
        pass
    return False


def get_aws_credentials():
    """Get AWS credentials secrets.

    Returns:
        Tuple[str, str, str]: AWS access key ID, AWS secret access key,
        and AWS default region.

    Raises:
        RuntimeError: If any of the required secrets are missing.
    """

    PROJECT_SECRETS_DIR = Path("./secrets")

    if _is_running_in_docker() or is_ci_environment():
        (
            aws_access_key_id,
            aws_secret_access_key,
            aws_default_region,
        ) = get_aws_credentials_from_docker_secrets()
    elif (
        os.getenv("AWS_ACCESS_KEY_ID")
        and os.getenv("AWS_SECRET_ACCESS_KEY")
        and os.getenv("AWS_DEFAULT_REGION")
    ):
        # Get AWS credentials from environment variables
        aws_access_key_id = os.getenv("AWS_ACCESS_KEY_ID")
        aws_secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")
        aws_default_region = os.getenv("AWS_DEFAULT_REGION")
    # check secrets folder exists with
    elif os.path.exists(PROJECT_SECRETS_DIR):
        (
            aws_access_key_id,
            aws_secret_access_key,
            aws_default_region,
        ) = _get_aws_credentials(PROJECT_SECRETS_DIR)

    else:
        logging.error("AWS credentials not found.")
        logging.error(
            "Please specify via environment variables, or localy in /secrets folder."
        )
        raise RuntimeError("AWS credentials not found.")

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
    if is_ci_environment():
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
