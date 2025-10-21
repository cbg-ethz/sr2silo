"""SILO schema fetching and loading functionality.

This module handles:
- Fetching SILO database configuration from the API
- Caching schema locally for performance and offline use
- Parsing SILO schema to extract field information
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import requests


class SiloSchema:
    """Manages SILO database schema fetching and caching."""

    def __init__(
        self, lapis_url: str, organism: str, cache_dir: Optional[Path] = None
    ) -> None:
        """Initialize SILO schema manager.

        Args:
            lapis_url: Base URL for the LAPIS API (e.g., 'https://lapis.wasap.genspectrum.org')
            organism: Organism identifier (e.g., 'covid', 'sc2')
            cache_dir: Optional directory for caching schema. If None, uses default.
        """
        self.lapis_url = lapis_url.rstrip("/")
        self.organism = organism
        self.cache_dir = cache_dir or Path.home() / ".sr2silo" / "cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._schema: Optional[Dict] = None

    @property
    def cache_file(self) -> Path:
        """Get the cache file path for this organism's schema."""
        return self.cache_dir / f"silo_schema_{self.organism}.json"

    def fetch_schema(self, use_cache: bool = True) -> Dict:
        """Fetch SILO database schema from API or cache.

        Args:
            use_cache: If True, try to load from cache first and save to cache on fetch.

        Returns:
            Dict containing the SILO database configuration.

        Raises:
            Exception: If fetching fails and no cache is available.
        """
        # Try cache first if enabled
        if use_cache and self.cache_file.exists():
            try:
                with self.cache_file.open("r") as f:
                    self._schema = json.load(f)
                logging.info(
                    f"Loaded SILO schema for {self.organism} from cache: {self.cache_file}"
                )
                return self._schema
            except (json.JSONDecodeError, IOError) as e:
                logging.warning(
                    f"Failed to load cached schema: {e}. Fetching from API..."
                )

        # Fetch from API
        url = f"{self.lapis_url}/info/database-config"
        headers = {"accept": "application/json"}

        try:
            logging.info(f"Fetching SILO schema from {url}")
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()
            self._schema = response.json()

            # Save to cache if enabled
            if use_cache:
                try:
                    with self.cache_file.open("w") as f:
                        json.dump(self._schema, f, indent=2)
                    logging.info(f"Cached SILO schema to {self.cache_file}")
                except IOError as e:
                    logging.warning(f"Failed to cache schema: {e}")

            return self._schema

        except (requests.exceptions.RequestException, Exception) as e:
            # Try cache as fallback
            if self.cache_file.exists():
                logging.warning(
                    f"API fetch failed: {e}. Attempting to use cached schema..."
                )
                try:
                    with self.cache_file.open("r") as f:
                        self._schema = json.load(f)
                    logging.info("Using cached schema as fallback")
                    return self._schema
                except (json.JSONDecodeError, IOError) as cache_error:
                    logging.error(f"Cache fallback also failed: {cache_error}")

            logging.error(f"Failed to fetch SILO schema: {e}")
            raise Exception(f"Failed to fetch SILO schema and no cache available: {e}")

    def get_metadata_fields(self) -> List[str]:
        """Extract list of metadata field names from schema.

        Returns:
            List of metadata field names defined in SILO schema.

        Raises:
            ValueError: If schema hasn't been fetched yet.
        """
        if self._schema is None:
            raise ValueError(
                "Schema not loaded. Call fetch_schema() first."
            )

        # Extract metadata fields from schema
        # SILO schema structure typically has a 'metadata' section
        metadata = self._schema.get("metadata", [])
        if isinstance(metadata, list):
            return [field.get("name") for field in metadata if "name" in field]
        elif isinstance(metadata, dict):
            return list(metadata.keys())
        else:
            logging.warning(f"Unexpected metadata format in SILO schema: {type(metadata)}")
            return []

    def get_nucleotide_sequences(self) -> List[str]:
        """Extract list of nucleotide sequence names from schema.

        Returns:
            List of nucleotide sequence segment names (e.g., ['main']).
        """
        if self._schema is None:
            raise ValueError("Schema not loaded. Call fetch_schema() first.")

        nucleotide_sequences = self._schema.get("nucleotideSequences", [])
        if isinstance(nucleotide_sequences, list):
            return [seq.get("name") for seq in nucleotide_sequences if "name" in seq]
        return []

    def get_genes(self) -> List[str]:
        """Extract list of gene names from schema.

        Returns:
            List of gene names (e.g., ['S', 'ORF1a', 'N']).
        """
        if self._schema is None:
            raise ValueError("Schema not loaded. Call fetch_schema() first.")

        genes = self._schema.get("genes", [])
        if isinstance(genes, list):
            return [gene.get("name") for gene in genes if "name" in gene]
        return []

    def clear_cache(self) -> None:
        """Remove cached schema file."""
        if self.cache_file.exists():
            self.cache_file.unlink()
            logging.info(f"Cleared cache file: {self.cache_file}")
