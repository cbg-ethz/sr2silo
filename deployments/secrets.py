#!/usr/bin/env python3
"""GPG credential management for sr2silo deployments.

Usage:
    # Load credentials into environment
    eval "$(python deployments/secrets.py)"

    # Or import in Python
    from deployments.secrets import load_credentials
    creds = load_credentials()
"""

import os
import subprocess
import sys
from pathlib import Path
from getpass import getpass


def get_passphrase() -> str:
    """Get GPG passphrase from environment or prompt."""
    if pp := os.environ.get("GPG_PASSPHRASE"):
        return pp
    return getpass("GPG passphrase: ")


def load_credentials() -> dict[str, str]:
    """Decrypt and return credentials as dict."""
    deploy_dir = Path(__file__).parent
    enc_file = deploy_dir / "secrets" / "credentials.enc"

    if not enc_file.exists():
        raise FileNotFoundError(f"No credentials at {enc_file}")

    passphrase = get_passphrase()

    result = subprocess.run(
        [
            "gpg",
            "--quiet",
            "--batch",
            "--passphrase-fd",
            "0",
            "--decrypt",
            str(enc_file),
        ],
        input=passphrase.encode(),
        capture_output=True,
    )

    if result.returncode != 0:
        raise ValueError(f"Decryption failed: {result.stderr.decode()}")

    creds = {}
    for line in result.stdout.decode().strip().split("\n"):
        if line and not line.startswith("#") and "=" in line:
            k, v = line.split("=", 1)
            creds[k.strip()] = v.strip()

    return creds


def load() -> None:
    """Load credentials into environment and print for shell eval."""
    creds = load_credentials()
    for k, v in creds.items():
        os.environ[k] = v
        print(f"export {k}={v}")
    print(f"# ✓ Loaded {len(creds)} credentials", file=sys.stderr)


if __name__ == "__main__":
    load()
