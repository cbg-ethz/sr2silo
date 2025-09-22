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

        # Fetch sequences
        print("Fetching sequences from WASAP GenSpectrum API...")
        data = client.get_sequences()

        if data:
            print("Request successful!")
            print(f"Response type: {type(data)}")

            # Pretty print the JSON response (first few items if it's a list)
            if isinstance(data, list) and len(data) > 0:
                print(f"\nFound {len(data)} sequences")
                print("First sequence:")
                print(json.dumps(data[0], indent=2))
            else:
                print("\nResponse data:")
                print(json.dumps(data, indent=2))
        else:
            print("Request failed!")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
