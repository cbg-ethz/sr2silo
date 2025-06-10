""" "Interactions with the Lapis API."""

from __future__ import annotations

import csv
import json
import logging
import uuid
from pathlib import Path
from typing import List, TypedDict

import requests

from sr2silo.config import is_ci_environment


class RequestUploadResponse(TypedDict):
    """Type definition for request-upload response from Lapis API."""

    fileId: str
    url: str


class SubmitResponse(TypedDict):
    """Type definition for submit response from Lapis API."""

    status: str
    message: str


class FileReference(TypedDict):
    """Type definition for file reference in fileMapping."""

    fileId: str
    name: str


class LapisClient:
    """Client for interacting with the Lapis API."""

    def __init__(self, token_url: str, submission_url: str, organism: str) -> None:
        """Initialize the Lapis client.

        Args:
            token_url: URL for authentication token endpoint
            submission_url: Base URL for submission endpoint
            organism: Organism identifier (e.g., 'sc2', 'sars-cov-2')
        """
        self.token_url = token_url
        self.submission_url = submission_url
        self.organism = organism
        self.is_ci_environment = is_ci_environment()
        self.token = None

    def authenticate(self, username: str, password: str) -> None:
        """Authenticate with the Lapis API."""

        if self.is_ci_environment is True:
            logging.info("CI environment detected. Using dummy token.")
            self.token = "dummy_token"
            return None

        response = requests.post(
            self.token_url,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data={
                "username": username,
                "password": password,
                "grant_type": "password",
                "client_id": "backend-client",
            },
        )

        if response.status_code == 200:
            self.token = response.json().get("access_token")
            return None
        else:
            raise Exception(
                f"Error: Unable to authenticate. Status code: {response.status_code},"
                f"Response: {response.text}"
            )

    def submit(
        self,
        group_id: int,
        metadata_file_path: Path,
        processed_file_path: Path | None = None,
        submission_id: str | None = None,
        data_use_terms_type: str = "OPEN",
    ) -> SubmitResponse:
        """Submit data to the Lapis API using pre-signed upload approach.

        Args:
            group_id: The group ID for the submission
            metadata_file_path: Path to the metadata TSV file
            processed_file_path: Path to the processed data file (e.g., .ndjson.zst)
            submission_id: Unique identifier for this submission (auto-generated if not provided)
            data_use_terms_type: Data use terms type (default: "OPEN")

        Returns:
            Submit response object

        Raises:
            Exception: If submission fails or authentication token is missing
        """
        if self.token is None:
            raise Exception(
                "Authentication required. Please call authenticate() first."
            )

        # Generate submission_id if not provided
        if submission_id is None:
            submission_id = str(uuid.uuid4())
            logging.info(f"Generated submission_id: {submission_id}")

        if self.is_ci_environment:
            logging.info("CI environment detected. Returning mock submit response.")
            return {
                "status": "success",
                "message": "Mock submission successful in CI environment",
            }

        # Step 1: Upload processed file via pre-signed URL if provided
        file_mapping_entries = {}

        if processed_file_path is not None:
            # Request pre-signed URL for processed file
            upload_responses = self.request_upload(group_id=group_id, numberFiles=1)

            if not upload_responses:
                raise Exception("Failed to get upload URLs for processed file")

            processed_upload = upload_responses[0]

            # Upload processed file to S3 using pre-signed URL
            try:
                with open(processed_file_path, "rb") as f:
                    upload_response = requests.put(
                        processed_upload["url"],
                        data=f,
                        headers={"Content-Type": "application/octet-stream"},
                    )
                    upload_response.raise_for_status()
                    logging.info(
                        f"Processed file uploaded successfully to S3: {processed_file_path.name}"
                    )

                # Add to file mapping
                file_mapping_entries["sequences"] = [
                    {
                        "fileId": processed_upload["fileId"],
                        "name": processed_file_path.name,
                    }
                ]
            except Exception as e:
                logging.error(f"Failed to upload processed file: {e}")
                raise Exception(f"Processed file upload failed: {e}")

        # Step 2: Create file mapping
        file_mapping = {submission_id: file_mapping_entries}

        # Generate a unique request ID
        request_id = str(uuid.uuid4())

        # Step 4: Submit to the API
        url = f"{self.submission_url}/{self.organism}/submit"
        params = {"groupId": group_id, "dataUseTermsType": data_use_terms_type}

        headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {self.token}",
            "x-request-id": request_id,
        }

        # Create empty sequence file for form data
        files = {
            "metadataFile": (
                "metadata.tsv",
                open(metadata_file_path, "rb"),
                "text/tab-separated-values",
            ),
            "sequenceFile": ("", "", "application/octet-stream"),  # Empty sequence file
            "fileMapping": ("", json.dumps(file_mapping), "application/json"),
        }

        response = None
        try:
            response = requests.post(url, params=params, headers=headers, files=files)
            response.raise_for_status()

            logging.info("Submission successful.")

            # Parse response - adjust based on actual API response structure
            try:
                response_data = response.json()
                return {
                    "status": "success",
                    "message": response_data.get(
                        "message", "Submission completed successfully"
                    ),
                }
            except ValueError:
                # If response is not JSON, create a basic response
                return {
                    "status": "success",
                    "message": "Submission completed successfully",
                }

        except requests.exceptions.HTTPError as e:
            logging.error(f"Error submitting data to Lapis: {e}")
            if response is not None:
                logging.error(f"Response: {response.text}")
                raise Exception(
                    f"Failed to submit: {response.status_code} - {response.text}"
                )
            else:
                raise Exception(f"Failed to submit: {e}")
        except requests.exceptions.RequestException as e:
            logging.error(f"Network error during submission: {e}")
            raise Exception(f"Network error: {e}")
        finally:
            # Close file handles
            for file_tuple in files.values():
                if hasattr(file_tuple[1], "close"):
                    file_tuple[1].close()

    def request_upload(
        self, group_id: int, numberFiles: int
    ) -> List[RequestUploadResponse]:
        """Request S3 pre-signed URLs to upload files.

        Args:
            group_id: The group ID for the upload request
            numberFiles: The number of files to request upload URLs for

        Returns:
            List of request-upload response objects containing fileId and pre-signed URL

        Raises:
            Exception: If the request fails or authentication token is missing
        """
        if self.token is None:
            raise Exception(
                "Authentication required. Please call authenticate() first."
            )

        if self.is_ci_environment:
            logging.info("CI environment detected. Returning mock upload response.")
            # Return mock response for CI environment
            return [
                {
                    "fileId": str(uuid.uuid4()).upper(),
                    "url": f"https://dummyendpoint.com/dummybucket/files/{uuid.uuid4()}?mock=true",
                }
                for _ in range(numberFiles)
            ]

        # Generate a unique request ID
        request_id = str(uuid.uuid4())

        # Construct the URL with query parameters
        url = "https://backend-wise-files.loculus.org/files/request-upload"
        params = {"groupId": group_id, "numberFiles": numberFiles}

        headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {self.token}",
            "x-request-id": request_id,
        }

        response = None
        try:
            response = requests.post(url, params=params, headers=headers, data="")
            response.raise_for_status()

            upload_responses = response.json()
            logging.info(f"Successfully requested {len(upload_responses)} upload URLs")

            return upload_responses

        except requests.exceptions.HTTPError as e:
            logging.error(f"Error requesting upload URLs: {e}")
            if response is not None:
                logging.error(f"Response: {response.text}")
                raise Exception(
                    f"Failed to request upload URLs: {response.status_code} - {response.text}"
                )
            else:
                raise Exception(f"Failed to request upload URLs: {e}")
        except requests.exceptions.RequestException as e:
            logging.error(f"Network error requesting upload URLs: {e}")
            raise Exception(f"Network error: {e}")


class Submission:
    """Submission-related utilities.
    Methods for generating placeholder FASTA files containing "NNN" sequences,
    and S3 links"""

    def __init__(self, fasta: str, s3_link: Path) -> None:
        """Initialize the Submission object."""
        self.fasta = fasta
        self.s3_link = s3_link

    @staticmethod
    def generate_placeholder_fasta(submission_ids: list[str]) -> str:
        """
        Generates a placeholder FASTA file for each submission ID with "NNN" as
        the sequence.
        """
        fasta_entries = []
        for submission_id in submission_ids:
            fasta_entries.append(f">{submission_id}")
            fasta_entries.append("NNN")  # Placeholder sequence
        return "\n".join(fasta_entries)

    @staticmethod
    def get_submission_ids_from_tsv(file_path: Path) -> list[str]:
        """
        Reads a TSV file and extracts submission IDs by parsing the
        "submissionId" column.
        """
        submission_ids = []
        with open(file_path, "r") as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")

            # Check if "submissionId" exists in the header
            if (
                reader.fieldnames is not None
                and "submissionId" not in reader.fieldnames
            ):
                raise ValueError(
                    'Error: "submissionId" column not found in the TSV file.'
                )

            # Extract submission IDs from the "submissionId" column
            for row in reader:
                submission_ids.append(row["submissionId"])

        return submission_ids

    @staticmethod
    def create_metadata_file(result_dir: Path) -> Path:
        """Create a metadata TSV file with the format containing only a date field.

        Args:
            result_dir: Directory where to save the metadata file

        Returns:
            Path to the created metadata file
        """
        result_dir_submission = result_dir / "submission"
        result_dir_submission.mkdir(parents=True, exist_ok=True)

        metadata_fp = result_dir_submission / "metadata.tsv"
        with metadata_fp.open("w") as f:
            # Write header with just the date field for now
            f.write("date\n")
            # Write empty data for now as requested
            f.write("\n")

        logging.info(f"Metadata file created at: {metadata_fp}")
        return metadata_fp
