""" "Interactions with the Loculus API."""

from __future__ import annotations

import csv
import json
import logging
import uuid
from datetime import date
from pathlib import Path
from typing import Dict, List, TypedDict

import requests
import zstandard as zstd

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
        processed_file_path: Path,
        nucleotide_alignment: Path,
        submission_id: str | None = None,
        data_use_terms_type: str = "OPEN",
    ) -> SubmitResponse:
        """Submit data to the Lapis API using pre-signed upload approach.

        Args:
            group_id: The group ID for the submission
            metadata_file_path: Path to the metadata TSV file
            processed_file_path: Path to the processed data file (e.g., .ndjson.zst)
            nucleotide_alignment: Path to nucleotide alignment file (.bam)
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

        # Step 1: Upload both files via pre-signed URLs
        file_mapping_entries = {}

        # Always upload both files
        files_to_upload = [
            ("processed", processed_file_path),
            ("nucleotide", nucleotide_alignment),
        ]

        # Request pre-signed URLs for both files
        upload_responses = self.request_upload(group_id=group_id, numberFiles=2)

        if len(upload_responses) != 2:
            raise Exception(
                f"Failed to get upload URLs for both files. Expected 2, got {len(upload_responses)}"
            )

        # Upload each file to S3 using pre-signed URLs
        for i, (file_type, file_path) in enumerate(files_to_upload):
            upload_info = upload_responses[i]

            try:
                # Get file size for Content-Length header
                file_size = file_path.stat().st_size
                logging.info(
                    f"Uploading {file_type} file: {file_path.name} "
                    f"({file_size} bytes)"
                )

                with open(file_path, "rb") as f:
                    upload_response = requests.put(
                        upload_info["url"],
                        data=f,
                        headers={
                            "Content-Type": "application/octet-stream",
                            "Content-Length": str(file_size),
                        },
                    )
                    upload_response.raise_for_status()
                    logging.info(
                        f"{file_type.capitalize()} file uploaded successfully to S3: "
                        f"{file_path.name} ({file_size} bytes)"
                    )

                # Add to file mapping with correct key names
                if file_type == "processed":
                    file_mapping_entries["siloReads"] = [
                        {
                            "fileId": upload_info["fileId"],
                            "name": file_path.name,
                        }
                    ]
                elif file_type == "nucleotide":
                    file_mapping_entries["nucleotideAlignment"] = [
                        {
                            "fileId": upload_info["fileId"],
                            "name": file_path.name,
                        }
                    ]

            except Exception as e:
                logging.error(f"Failed to upload {file_type} file: {e}")
                raise Exception(f"{file_type.capitalize()} file upload failed: {e}")

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
                metadata_file_path.name,
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
    def create_metadata_file(
        processed_file: Path, count_reads: bool = False
    ) -> tuple[Path, str]:
        """Create a metadata TSV file with the required submissionId header.

            The metadata file will be saved in a 'submission' subdirectory of
            its parent directory.

        Args:
            processed_file: Path to the processed file.
            count_reads: Whether to include a countReads column (default: False)

        Returns:
            Tuple of (Path to the created metadata file, submission ID)
        """
        # Ensure the submission directory exists
        result_dir = processed_file.parent
        result_dir_submission = result_dir / "submission"
        result_dir_submission.mkdir(parents=True, exist_ok=True)

        # Generate a submission id and current date
        submission_id = str(uuid.uuid4())
        submission_date = date.today().isoformat()

        # Count reads if requested
        count = None
        if count_reads:
            count = Submission.count_reads(processed_file)
            logging.info(f"Counted {count} reads in {processed_file}")
        else:
            logging.debug("Skipping read count as countReads is False")

        # extract metadata from the processed file
        metadata = Submission.parse_metadata(processed_file)
        logging.debug(f"Extracted metadata: {metadata}")

        # Create the metadata file path

        metadata_fp = result_dir_submission / f"metadata_{submission_id}.tsv"
        with metadata_fp.open("w") as f:
            # Prepare all column names
            columns = ["submissionId", "date"]
            values = [submission_id, submission_date]

            # Add metadata fields as columns
            for key, value in metadata.items():
                columns.append(key)
                values.append(str(value) if value is not None else "")

            # Add countReads column if requested
            if count is not None:
                columns.append("countReads")
                values.append(str(count))

            # Write header row
            f.write("\t".join(columns) + "\n")
            # Write data row
            f.write("\t".join(values) + "\n")
        logging.info(
            f"Metadata file created at: {metadata_fp} "
            f"with submission ID: {submission_id}"
        )
        return metadata_fp, submission_id

    @staticmethod
    def parse_metadata(silo_input: Path) -> Dict:
        """Parses the metadata from a silo input .ndjson.zstd or .ndjson
        returning all metadata fields but readId as a dictoriary with keys
        in camel case.

        Assumptions:
         - the metaddata is stored in the root of the object under the keys
         - each read has the same metadata, up to readId
        """
        # check if the file is compressed or not
        if silo_input.suffix == ".zst" or silo_input.name.endswith(".zstd"):
            # open the file accordingly (compressed)
            with open(silo_input, "rb") as f:
                dctx = zstd.ZstdDecompressor()
                with dctx.stream_reader(f) as reader:
                    # read the first line by reading chunks until we find a newline
                    buffer = b""
                    chunk_size = 1024
                    first_line = ""
                    while True:
                        chunk = reader.read(chunk_size)
                        if not chunk:
                            # If no more data and no newline found,
                            # use the entire buffer
                            first_line = buffer.decode("utf-8").strip()
                            break
                        buffer += chunk
                        if b"\n" in buffer:
                            first_line = (
                                buffer.split(b"\n", 1)[0].decode("utf-8").strip()
                            )
                            break
        else:
            # open the file accordingly (uncompressed)
            with open(silo_input, "r") as f:
                # read the first line
                first_line = f.readline().strip()

        # parse the JSON and extract the metadata - look for fields sample_id, batch_id
        # location_code, location_name, sampling_date, sr2silo_version
        data = json.loads(first_line)

        # Extract specific metadata fields
        metadata_fields = [
            "sample_id",
            "batch_id",
            "location_code",
            "location_name",
            "sampling_date",
            "sr2silo_version",
        ]

        metadata = {}
        for field in metadata_fields:
            if field in data:
                metadata[field] = data[field]
            else:
                logging.warning(
                    f"Metadata field '{field}' not found in the input data."
                )
                metadata[field] = None

        # convert keys to camel case
        def to_camel_case(snake_str):
            """Convert snake_case string to camelCase string."""
            components = snake_str.split("_")
            return components[0] + "".join(word.capitalize() for word in components[1:])

        camel_case_metadata = {}
        for key, value in metadata.items():
            if key != "read_id":  # exclude readId as specified
                camel_case_key = to_camel_case(key)
                camel_case_metadata[camel_case_key] = value

        # return the metadata dictionary
        return camel_case_metadata

    @staticmethod
    def count_reads(silo_input: Path) -> int:
        """Counts the number of reads in a silo input .ndjson.zstd or .ndjson file.

        Assumption: each line in the file corresponds to one read.
        """

        if silo_input.suffix == ".zst" or silo_input.name.endswith(".zstd"):
            # Handle compressed files
            with open(silo_input, "rb") as f:
                dctx = zstd.ZstdDecompressor()
                with dctx.stream_reader(f) as reader:
                    count = 0
                    chunk_size = 8192  # Read in 8KB chunks
                    last_byte = None
                    while True:
                        chunk = reader.read(chunk_size)
                        if not chunk:
                            break
                        count += chunk.count(b"\n")
                        last_byte = chunk[-1] if len(chunk) > 0 else last_byte
            # Only add 1 if the last byte is not a newline
            if last_byte is not None and last_byte != ord(b"\n"):
                return count + 1
            else:
                return count
        else:
            # Handle uncompressed files
            count = 0
            last_byte = None
            with open(silo_input, "rb") as f:
                chunk_size = 8192  # Read in 8KB chunks
                while True:
                    chunk = f.read(chunk_size)
                    if not chunk:
                        break
                    count += chunk.count(b"\n")
                    last_byte = chunk[-1] if len(chunk) > 0 else last_byte
            # Only add 1 if the last byte is not a newline
            if last_byte is not None and last_byte != ord(b"\n"):
                return count + 1
            else:
                return count
