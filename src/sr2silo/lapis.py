""""Interactions with the Lapis API."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import requests

# TODO: move to environment variables
KEYCLOAK_TOKEN_URL = (
    "https:"
    "//authentication-wise-seqs.loculus.org"
    "/realms/loculus/protocol/openid-connect/token"
)
SUBMISSION_URL = (
    "https:"
    "//backend-wise-seqs.loculus.org"
    "/test/submit?groupId={group_id}&dataUseTermsType=OPEN"
)


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


def get_submission_ids_from_tsv(file_path: str) -> list[str]:
    """
    Reads a TSV file and extracts submission IDs by parsing the "submissionId" column.
    """
    submission_ids = []
    with open(file_path, "r") as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter="\t")

        # Check if "submissionId" exists in the header
        if reader.fieldnames is not None and "submissionId" not in reader.fieldnames:
            raise ValueError('Error: "submissionId" column not found in the TSV file.')

        # Extract submission IDs from the "submissionId" column
        for row in reader:
            submission_ids.append(row["submissionId"])

    return submission_ids


def get_loculus_authentication_token(username: str, password: str) -> str:
    """
    Sends a request to the Keycloak authentication server to obtain a token.
    """
    response = requests.post(
        KEYCLOAK_TOKEN_URL,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        data={
            "username": username,
            "password": password,
            "grant_type": "password",
            "client_id": "backend-client",
        },
    )

    if response.status_code == 200:
        return response.json().get("access_token")
    else:
        raise Exception(
            f"Error: Unable to authenticate. Status code: {response.status_code},"
            f"Response: {response.text}"
        )


def _submit(
    authentication_token: str, group_id: int, tsv_path: str, fasta_path: str
) -> None:
    """
    Submits the metadata and sequence files to Loculus via a POST request.
    """
    submission_url = SUBMISSION_URL.format(group_id=group_id)

    with open(tsv_path, "rb") as tsv_file, open(fasta_path, "rb") as fasta_file:
        response = requests.post(
            submission_url,
            headers={
                "Authorization": f"Bearer {authentication_token}",
                "accept": "application/json",
            },
            files={"metadataFile": tsv_file, "sequenceFile": fasta_file},
        )

    if response.status_code == 200:
        print("Upload successful.")
        print(
            "You can approve the upload for release at:\n\n"
            "https://wise-seqs.loculus.org/salmonella/submission/1/review"
        )
    else:
        raise Exception(
            f"Error: Unable to submit. Status code: {response.status_code}, "
            f"Response: {response.text}"
        )


def submit(input_fp: Path, username: str, password: str, group_id: int) -> None:
    """
    Upload the a metadata tsv file to a loculus instance.
    """

    authentication_token = get_loculus_authentication_token(username, password)
    submission_ids = get_submission_ids_from_tsv(str(input_fp))
    placeholder_fasta_str = generate_placeholder_fasta(submission_ids)

    # Write the placeholder FASTA to a temporary file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as fasta_file:
        fasta_file.write(placeholder_fasta_str.encode("utf-8"))
        placeholder_tmp_path = fasta_file.name

    _submit(authentication_token, group_id, str(input_fp), placeholder_tmp_path)
