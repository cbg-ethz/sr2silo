""" "Interactions with the Loculus API."""

from __future__ import annotations

import csv
import json
import logging
import uuid
from datetime import date
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


class LoculusClient:
    """Client for interacting with the Loculus API."""

    def __init__(self, token_url: str, submission_url: str, organism: str) -> None:
        """Initialize the Loculus client.

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
        """Authenticate with the Loculus API."""

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
            submission_id: Unique identifier for this submission
                           (auto-generated if not provided)
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

        # Use provided submission_id or generate if not provided
        if submission_id is None:
            submission_id = str(uuid.uuid4())
            logging.info(f"Generated submission_id: {submission_id}")
        else:
            logging.info(f"Using provided submission_id: {submission_id}")

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
                # Get file size for Content-Length header
                file_size = processed_file_path.stat().st_size
                logging.info(
                    f"Uploading processed file: {processed_file_path.name} "
                    f"({file_size} bytes)"
                )

                with open(processed_file_path, "rb") as f:
                    upload_response = requests.put(
                        processed_upload["url"],
                        data=f,
                        headers={
                            "Content-Type": "application/octet-stream",
                            "Content-Length": str(file_size),
                        },
                    )
                    upload_response.raise_for_status()
                    logging.info(
                        f"Processed file uploaded successfully to S3: "
                        f"{processed_file_path.name} ({file_size} bytes)"
                    )

                # Add to file mapping
                file_mapping_entries["silo_reads"] = [
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

        # Create files for form data - no sequence file for sc2 organism
        files = {
            "metadataFile": (
                "metadata.tsv",
                open(metadata_file_path, "rb"),
                "text/tab-separated-values",
            ),
            "fileMapping": ("", json.dumps(file_mapping), "application/json"),
        }

        response = None
        try:
            response = requests.post(url, params=params, headers=headers, files=files)
            response.raise_for_status()

            logging.info("Submission successful.")
            logging.debug(f"Response status: {response.status_code}")
            logging.debug(f"Response headers: {response.headers}")
            logging.debug(f"Response text: {response.text}")

            # Parse response
            response_data = response.json()
            accession = response_data[0].get("accession", "N/A")
            submission_id = response_data[0].get("submissionId", "N/A")
            return {
                "status": "success",
                "message": f"Submission completed successfully."
                f"Accession: {accession}, Submission ID: {submission_id}",
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
                    "url": f"""https://
                    dummyendpoint.com/dummybucket/files/{uuid.uuid4()}?mock=true""",
                }
                for _ in range(numberFiles)
            ]

        # Generate a unique request ID
        request_id = str(uuid.uuid4())

        # Construct the URL with query parameters
        url = f"{self.submission_url}/files/request-upload"
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
                    f"""Failed to request upload URLs:
                    {response.status_code} - {response.text}"""
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
    def create_metadata_file(result_dir: Path) -> tuple[Path, str]:
        """Create a metadata TSV file with the required submissionId header.

        Args:
            result_dir: Directory where to save the metadata file

        Returns:
            Tuple of (Path to the created metadata file, submission ID)
        """
        result_dir_submission = result_dir / "submission"
        result_dir_submission.mkdir(parents=True, exist_ok=True)

        submission_id = str(uuid.uuid4())
        metadata_fp = result_dir_submission / "metadata.tsv"
        with metadata_fp.open("w") as f:
            # Write header with required submissionId field and optional date field
            f.write("submissionId\tdate\n")
            # Write a sample entry with the generated UUID for submissionId
            today = date.today().isoformat()
            f.write(f"{submission_id}\t{today}\n")

        logging.info(
            f"Metadata file created at: {metadata_fp} "
            f"with submission ID: {submission_id}"
        )
        return metadata_fp, submission_id


def refererenceGenome():
    """Fetch reference genome from the Lapis `sample/referenceGenome` endpoint"""


def _reference_json_to_fasta(
    reference_json_string: str, nucleotide_out_fp: Path, amino_acid_out_fp: Path
) -> None:
    """Convert a reference JSON from `sample/referenceGenome` endpoint
      to separate nucleotide and amino acid reference FASTA files.

    Args:
        reference_json_string: JSON string containing reference sequences with
                              'nucleotideSequences' and 'genes' sections
        nucleotide_output_path: Path to the output nucleotide FASTA file
        amino_acid_output_path: Path to the output amino acid FASTA file

    Returns:
        None
    """
    # Parse the JSON string
    reference_data = json.loads(reference_json_string)

    # Create nucleotide FASTA file
    with nucleotide_out_fp.open("w") as nuc_file:
        for nuc_seq in reference_data.get("nucleotideSequences", []):
            seq_name = nuc_seq.get("name", "main")
            sequence = nuc_seq.get("sequence", "")
            nuc_file.write(f">{seq_name}\n{sequence}\n")

    logging.info(f"Nucleotide FASTA file created at: {nucleotide_out_fp}")

    # Create amino acid FASTA file
    with amino_acid_out_fp.open("w") as aa_file:
        for gene in reference_data.get("genes", []):
            gene_name = gene.get("name", "")
            sequence = gene.get("sequence", "")
            # Remove stop codon asterisk if present
            sequence = sequence.rstrip("*")
            aa_file.write(f">{gene_name}\n{sequence}\n")

    logging.info(f"Amino acid FASTA file created at: {amino_acid_out_fp}")
