"""Functionality to interact with AWS S3."""

from __future__ import annotations

import os

import boto3

# Access AWS credentials from environment variables
AWS_ACCESS_KEY_ID = os.getenv("AWS_ACCESS_KEY_ID")
AWS_SECRET_ACCESS_KEY = os.getenv("AWS_SECRET_ACCESS_KEY")
AWS_DEFAULT_REGION = os.getenv("AWS_DEFAULT_REGION")


# Create an S3 client
s3 = boto3.client(
    "s3",
    aws_access_key_id=AWS_ACCESS_KEY_ID,
    aws_secret_access_key=AWS_SECRET_ACCESS_KEY,
    region_name=AWS_DEFAULT_REGION,
)

# bucket_name = 'vpipe-output'
# kp3_mutations = 'mut_def/kp23.yaml'
# xec_mutations = "mut_def/xec.yaml"


def upload_file_to_s3(file_name, bucket_name, object_name=None):
    """Upload a file to an S3 bucket.

    Args:
        file_name (str): Path to the file to upload.
        bucket_name (str): Bucket to upload to.
        object_name (str, optional): S3 object name.
                        If not specified then file_name is used.
    """
    if object_name is None:
        object_name = os.path.basename(file_name)

    try:
        s3.upload_file(file_name, bucket_name, object_name)
        print(f"File {file_name} uploaded to {bucket_name}/{object_name}")
    except Exception as e:
        print(f"Error uploading file to S3: {e}")


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

    try:
        s3.download_file(bucket_name, object_name, file_name)
        print(f"File {object_name} downloaded from {bucket_name} to {file_name}")
    except Exception as e:
        print(f"Error downloading file from S3: {e}")


def list_files_in_bucket(bucket_name):
    """List all files in an S3 bucket.

    Args:
        bucket_name (str): The name of the S3 bucket.
    """
    try:
        response = s3.list_objects_v2(Bucket=bucket_name)
        if "Contents" in response:
            files = [obj["Key"] for obj in response["Contents"]]
            return files
        else:
            return []
    except Exception as e:
        print(f"Error listing files in bucket {bucket_name}: {e}")
        return []
