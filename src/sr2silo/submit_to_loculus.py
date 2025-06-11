"""Submit processd data to Loculus Silo and upload to S3."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from sr2silo.config import (

    get_frontend_url,
    get_keycloak_token_url,
    get_mock_urls,
    get_organism,
    get_submission_url,
    is_ci_environment,
)
from sr2silo.silo import LapisClient, Submission
from sr2silo.storage import upload_file_to_s3


def load_config(config_file: Path) -> dict:
    """Load a JSON configuration file."""
    try:
        with config_file.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"Config file not found: {config_file}")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from config file: {config_file} - {e}")
        raise


def make_submission_file(result_dir: Path, srLink: str) -> Path:
    """Create a submission file with the given S3 link.

    Args:
        result_dir (Path): The directory to save the submission file.
        srLink (str): The S3 link to include in the submission file.

    Returns:
        Path: The path to the created submission file.
    """
    result_dir_submission = result_dir / "submission"
    result_dir_submission.mkdir(parents=True, exist_ok=True)

    submission_metadata_fp = result_dir_submission / "metadata.tsv"
    with submission_metadata_fp.open("w") as f:
        f.write("submissionId\ts3Link\tversionComment\n")
        f.write(f"001\t{srLink}\t\n")
    logging.info(f"Submission metadata saved to: {submission_metadata_fp}")

    return submission_metadata_fp


def upload_to_s3(aligned_reads_fp: Path, sample_id: str) -> str:
    """Upload a file to S3 bucket.

    Args:
        aligned_reads_fp (Path): The file to upload.
        sample_id (str): Sample ID used in the S3 filename.

    Returns:
        str: The S3 link to the uploaded file.
    """
    logging.info(f"Uploading to S3: {aligned_reads_fp}")
    suffix = aligned_reads_fp.suffix
    s3_file_name = f"{sample_id}.{suffix}"
    s3_bucket = "sr2silo01"  # TODO : Make this configurable
    s3_link = f"s3://{s3_bucket}/{s3_file_name}"
    upload_file_to_s3(aligned_reads_fp, s3_bucket, s3_file_name)
    return s3_link



def submit_to_silo(result_dir: Path, processed_file: Path) -> bool:
    """Submit data to SILO using the new pre-signed upload approach.

    Args:
        result_dir (Path): Directory where to save submission files.
        processed_file (Path): Path to the processed .ndjson.zst file to upload.

    Returns:
        bool: True if submission was successful, False otherwise.
    """
    logging.info("Submitting to Loculus...")

    # Create the new metadata file format and get the submission ID
    metadata_fp, submission_id = Submission.create_metadata_file(result_dir)

    if is_ci_environment():
        logging.info(
            "CI environment active; using mock URLs but executing submission mechanics."
        )
        # Get mock URLs for CI environment, then continue with LAPIS submission
        KEYCLOAK_TOKEN_URL, SUBMISSION_URL = get_mock_urls()
    else:
        # Get URLs from environment or use defaults
        KEYCLOAK_TOKEN_URL = get_keycloak_token_url()
        SUBMISSION_URL = get_submission_url()
        logging.info(f"Using Keycloak URL: {KEYCLOAK_TOKEN_URL}")
        logging.info(f"Using submission URL: {SUBMISSION_URL}")

    # Get organism configuration
    organism = get_organism()
    logging.info(f"Using organism: {organism}")

    try:
        # Create client with organism parameter
        client = LapisClient(KEYCLOAK_TOKEN_URL, SUBMISSION_URL, organism)
        client.authenticate(username="testuser", password="testuser")

        # Submit using new API with both metadata and processed file
        response = client.submit(
            group_id=1,
            metadata_file_path=metadata_fp,
            processed_file_path=processed_file,
            submission_id=submission_id,
        )

        if response["status"] == "success":
            logging.info(response["message"])

            # Get the frontend URL from config
            frontend_url = get_frontend_url()
            logging.info(
                f"You can approve the upload for release at: "
                f"{frontend_url}/{organism}/submission/1/review"
            )
            return True
        else:
            logging.error(f"Submission failed: {response}")
            return False

    except Exception as e:
        logging.error(f"Error during submission: {e}")
        return False
