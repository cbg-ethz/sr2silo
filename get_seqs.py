from __future__ import annotations

import json
import os
import uuid

import requests


class GenSpectrumClient:
    """Client for interacting with the WASAP GenSpectrum API."""

    def __init__(
        self, token_url: str, api_base_url: str, organism: str = "covid"
    ) -> None:
        """Initialize the GenSpectrum client.

        Args:
            token_url: URL for authentication token endpoint
            api_base_url: Base URL for the API
            organism: Organism identifier (default: "covid")
        """
        self.token_url = token_url
        self.api_base_url = api_base_url
        self.organism = organism
        self.token = None

    def authenticate(self, username: str, password: str) -> None:
        """Authenticate with the API using username/password."""
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
            print("Authentication successful!")
        else:
            raise Exception(
                f"Authentication failed. Status code: {response.status_code}, "
                f"Response: {response.text}"
            )

    def get_sequences(self, **params):
        """
        Fetch sequences from the GenSpectrum API

        Args:
            **params: Additional query parameters for the API request
        """
        if self.token is None:
            raise Exception(
                "Authentication required. Please call authenticate() first."
            )

        url = f"{self.api_base_url}/backend/{self.organism}/get-sequences"

        # Generate a unique request ID
        request_id = str(uuid.uuid4())

        headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {self.token}",
            "x-request-id": request_id,
        }

        try:
            response = requests.get(url, headers=headers, params=params)

            print(f"Status code: {response.status_code}")
            print(f"Request ID: {request_id}")

            response.raise_for_status()

            # Return the JSON response
            return response.json()

        except requests.exceptions.RequestException as e:
            print(f"Error making request: {e}")
            if hasattr(e, "response") and e.response is not None:
                print(f"Response text: {e.response.text}")
            return None

    def get_original_metadata(
        self, statuses_filter: str = "APPROVED_FOR_RELEASE", **params
    ):
        """
        Fetch original metadata from the GenSpectrum API

        Args:
            statuses_filter: Filter by status (default: "APPROVED_FOR_RELEASE")
            **params: Additional query parameters for the API request
        """
        if self.token is None:
            raise Exception(
                "Authentication required. Please call authenticate() first."
            )

        url = f"{self.api_base_url}/backend/{self.organism}/get-original-metadata"

        # Generate a unique request ID
        request_id = str(uuid.uuid4())

        headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {self.token}",
            "x-request-id": request_id,
        }

        # Add the statusesFilter parameter
        params["statusesFilter"] = statuses_filter

        try:
            response = requests.get(url, headers=headers, params=params)

            print(f"Status code: {response.status_code}")
            print(f"Request ID: {request_id}")
            print(f"URL: {response.url}")
            print(
                f"Content-Type: {response.headers.get('Content-Type', 'Not specified')}"
            )
            print(
                f"Content-Length: {response.headers.get('Content-Length', 'Not specified')}"
            )

            response.raise_for_status()

            # Debug: Print raw response first
            raw_text = response.text
            print(f"Raw response (first 500 chars): {raw_text[:500]}")
            print(f"Raw response (last 500 chars): {raw_text[-500:]}")

            # Try to parse as JSON
            try:
                return response.json()
            except json.JSONDecodeError as json_error:
                print(f"JSON decode error: {json_error}")
                print("Response might be NDJSON or multiple JSON objects")

                # Try to parse as NDJSON (one JSON object per line)
                lines = raw_text.strip().split("\n")
                print(f"Found {len(lines)} lines in response")

                if len(lines) > 0:
                    try:
                        # Try to parse first line as JSON
                        first_obj = json.loads(lines[0])
                        print(f"First line parsed successfully: {type(first_obj)}")
                        if isinstance(first_obj, dict):
                            print(f"Keys in first object: {list(first_obj.keys())}")

                        # Return all parsed objects
                        parsed_objects = []
                        for i, line in enumerate(lines):
                            if line.strip():
                                try:
                                    parsed_objects.append(json.loads(line))
                                except json.JSONDecodeError:
                                    print(f"Failed to parse line {i}: {line[:100]}...")
                                    break

                        print(f"Successfully parsed {len(parsed_objects)} JSON objects")
                        return parsed_objects

                    except json.JSONDecodeError:
                        print("Could not parse as NDJSON either")
                        return raw_text
                else:
                    return raw_text

        except requests.exceptions.RequestException as e:
            print(f"Error making request: {e}")
            if hasattr(e, "response") and e.response is not None:
                print(f"Response text: {e.response.text}")
            return None


def main():
    """
    Main function to execute the API request
    """
    # Configuration from your config.yaml
    token_url = "https://auth.db.wasap.genspectrum.org/realms/loculus/protocol/openid-connect/token"
    api_base_url = "https://api.db.wasap.genspectrum.org"
    organism = "covid"

    # Get credentials from environment variables (more secure)
    username = os.getenv("LOCULUS_USERNAME")
    password = os.getenv("LOCULUS_PASSWORD")

    if not username or not password:
        print(
            "Error: Please set LOCULUS_USERNAME and LOCULUS_PASSWORD environment variables"
        )
        print("Example:")
        print("  export LOCULUS_USERNAME='your_username'")
        print("  export LOCULUS_PASSWORD='your_password'")
        return

    # Initialize client
    client = GenSpectrumClient(token_url, api_base_url, organism)

    try:
        # Authenticate
        print("Authenticating...")
        client.authenticate(username, password)

        print("=" * 80)
        print("EXPLORING: get-original-metadata endpoint")
        print("=" * 80)

        # Fetch original metadata
        print("Fetching original metadata from WASAP GenSpectrum API...")
        metadata = client.get_original_metadata()

        # set of sample_ids
        sample_ids = set()

        if metadata:
            print("Original metadata request successful!")
            print(f"Response type: {type(metadata)}")

            # Pretty print the JSON response (first few items if it's a list)
            if isinstance(metadata, list) and len(metadata) > 0:
                print(f"\nFound {len(metadata)} metadata entries")
                print("\n" + "=" * 80)
                print("COMPLETE LIST OF ALL METADATA ENTRIES:")
                print("=" * 80)

                for i, entry in enumerate(metadata, 1):
                    print(f"\n--- Entry {i:3d} ---")
                    print(f"Accession: {entry.get('accession', 'N/A')}")
                    print(f"Version: {entry.get('version', 'N/A')}")
                    print(f"Submitter: {entry.get('submitter', 'N/A')}")
                    print(f"Is Revocation: {entry.get('isRevocation', 'N/A')}")

                    if entry.get("isRevocation"):
                        print("REVOCATION ENTRY - No original metadata available")
                    else:
                        orig_meta = entry.get("originalMetadata", {})
                        sample_id = orig_meta.get("sampleId")
                        if sample_id:
                            sample_ids.add(sample_id)
                        print(f"Sample ID: {orig_meta.get('sampleId', 'N/A')}")
                        print(f"Batch ID: {orig_meta.get('batchId', 'N/A')}")
                        print(
                            f"Location: {orig_meta.get('locationName', 'N/A')} ({orig_meta.get('locationCode', 'N/A')})"
                        )
                        print(f"Sampling Date: {orig_meta.get('samplingDate', 'N/A')}")
                        print(f"Submission Date: {orig_meta.get('date', 'N/A')}")
                        print(f"Read Count: {orig_meta.get('countSiloReads', 'N/A')}")
                        print(
                            f"sr2silo Version: {orig_meta.get('sr2siloVersion', 'N/A')}"
                        )

                print("\n" + "=" * 80)
                print(f"TOTAL: {len(metadata)} entries displayed")

                print(f"UNIQUE SAMPLE IDs found: {len(sample_ids)}")
                print("=" * 80)
            else:
                print("\nMetadata response:")
                print(json.dumps(metadata, indent=2))
        else:
            print("Original metadata request failed!")

        print("#     print('Request failed!')")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
