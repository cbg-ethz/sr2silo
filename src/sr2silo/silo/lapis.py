""""Interactions with the Lapis API."""

from __future__ import annotations

import csv
import logging
import tempfile
from pathlib import Path

import requests

from sr2silo.config import is_ci_environment


class LapisClient:
    """Client for interacting with the Lapis API."""

    def __init__(self, token_url: str, submission_url: str) -> None:
        """Initialize the Lapis client."""
        self.token_url = token_url
        self.submission_url = submission_url
        self.is_ci_environment = is_ci_environment
        self.token = None

    def authenticate(self, username: str, password: str) -> None:
        """Authenticate with the Lapis API."""

        if self.is_ci_environment is True:
            logging.info("CI environment detected. Using dummy token.")
            self.token = "dummy_token"

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

    def submit(self, group_id: int, data: dict) -> requests.Response:
        """Submit data to the Lapis API."""

        if self.is_ci_environment is True:
            logging.info("Running in CI environment, skipping actual submission.")
            return requests.Response()

        url = self.submission_url.format(group_id=group_id)

        # Write the placeholder FASTA to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as fasta_file:
            fasta_file.write(data["fasta"].encode("utf-8"))
            placeholder_tmp_path = fasta_file.name

        with open(data["input_fp"], "rb") as tsv_file, open(
            placeholder_tmp_path, "rb"
        ) as fasta_file:
            response = requests.post(
                url,
                headers={
                    "Authorization": f"Bearer {self.token}",
                    "accept": "application/json",
                },
                files={"metadataFile": tsv_file, "sequenceFile": fasta_file},
            )

        response.raise_for_status()

        if response.status_code == 200:
            logging.info("Upload successful.")
            logging.info(
                "You can approve the upload for release at:\n\n"
                "https://wise-seqs.loculus.org/salmonella/submission/1/review"
            )
        else:
            error_message = (
                f"Error: Unable to submit. Status code: {response.status_code}, "
                f"Response: {response.text}"
            )
            logging.error(error_message)
            raise Exception(error_message)
        return response


class Submission:
    """Submission-related utilities.
    Methods for generating placeholder FASTA files containing "NNN" sequences,
    and S3 links"""

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
