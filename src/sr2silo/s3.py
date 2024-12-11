"""Functionality to interact with AWS S3."""

from __future__ import annotations

import os
import shutil

import boto3
from botocore.exceptions import NoCredentialsError


def compress_file(input_file, output_file):
    """Compress a file using BZ2 compression."""
    with open(input_file, "rb") as f_in:
        with open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def get_s3_client():
    """Get an S3 client using AWS credentials from environment variables."""

    # Get AWS credentials from environment variables
    aws_access_key_id = os.getenv("AWS_ACCESS_KEY_ID")
    aws_secret_access_key = os.getenv("AWS_SECRET_ACCESS_KEY")
    aws_default_region = os.getenv("AWS_DEFAULT_REGION")

    # Create an S3 client
    s3_client = boto3.client(
        "s3",
        aws_access_key_id=aws_access_key_id,
        aws_secret_access_key=aws_secret_access_key,
        region_name=aws_default_region,
    )

    return s3_client


def upload_file_to_s3(file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket"""
    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = file_name

    # get the s3 client
    s3_client = get_s3_client()

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
