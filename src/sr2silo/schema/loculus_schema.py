"""Loculus schema loading and management.

This module handles:
- Loading Loculus configuration from local YAML or dict
- Parsing Loculus schema to extract field information
- Managing organism-specific metadata field definitions
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import yaml
except ImportError:
    yaml = None


class LoculusSchema:
    """Manages Loculus schema loading and field extraction."""

    def __init__(self, organism: str, config: Optional[Dict[str, Any]] = None) -> None:
        """Initialize Loculus schema manager.

        Args:
            organism: Organism identifier (e.g., 'covid', 'sc2')
            config: Optional pre-loaded configuration dict. If None, must call load_from_file.
        """
        self.organism = organism
        self._config = config
        self._metadata_fields: Optional[List[Dict[str, Any]]] = None

    @classmethod
    def from_file(cls, organism: str, config_file: Path) -> "LoculusSchema":
        """Create LoculusSchema from a YAML configuration file.

        Args:
            organism: Organism identifier
            config_file: Path to YAML configuration file

        Returns:
            LoculusSchema instance with loaded configuration

        Raises:
            ImportError: If PyYAML is not installed
            FileNotFoundError: If config file doesn't exist
            ValueError: If organism not found in config
        """
        if yaml is None:
            raise ImportError("PyYAML is required to load YAML files. Install with: pip install pyyaml")

        if not config_file.exists():
            raise FileNotFoundError(f"Config file not found: {config_file}")

        with config_file.open("r") as f:
            full_config = yaml.safe_load(f)

        # Extract organism-specific config
        organisms_config = full_config.get("organisms", {})
        if organism not in organisms_config:
            raise ValueError(
                f"Organism '{organism}' not found in config. Available: {list(organisms_config.keys())}"
            )

        organism_config = organisms_config[organism]
        return cls(organism=organism, config=organism_config)

    @classmethod
    def from_dict(cls, organism: str, config: Dict[str, Any]) -> "LoculusSchema":
        """Create LoculusSchema from a dictionary configuration.

        Args:
            organism: Organism identifier
            config: Configuration dictionary (organism-specific portion)

        Returns:
            LoculusSchema instance with loaded configuration
        """
        return cls(organism=organism, config=config)

    def get_metadata_fields(self) -> List[str]:
        """Extract list of metadata field names from Loculus schema.

        Returns:
            List of metadata field names defined in Loculus schema.

        Raises:
            ValueError: If configuration hasn't been loaded.
        """
        if self._config is None:
            raise ValueError("Configuration not loaded.")

        # Parse metadata fields from schema
        schema = self._config.get("schema", {})
        metadata = schema.get("metadata", [])

        if not isinstance(metadata, list):
            logging.warning(f"Expected list for metadata, got {type(metadata)}")
            return []

        return [field.get("name") for field in metadata if isinstance(field, dict) and "name" in field]

    def get_metadata_field_info(self) -> List[Dict[str, Any]]:
        """Get detailed information about metadata fields.

        Returns:
            List of metadata field definitions with name, type, displayName, etc.
        """
        if self._config is None:
            raise ValueError("Configuration not loaded.")

        schema = self._config.get("schema", {})
        metadata = schema.get("metadata", [])

        if not isinstance(metadata, list):
            return []

        return [field for field in metadata if isinstance(field, dict)]

    def get_field_type(self, field_name: str) -> Optional[str]:
        """Get the type of a specific metadata field.

        Args:
            field_name: Name of the field

        Returns:
            Field type (e.g., 'string', 'date') or None if not found
        """
        field_info = self.get_metadata_field_info()
        for field in field_info:
            if field.get("name") == field_name:
                return field.get("type")
        return None

    def get_required_fields(self) -> List[str]:
        """Get list of required metadata fields.

        Returns:
            List of required field names
        """
        field_info = self.get_metadata_field_info()
        # Fields without explicit 'required: false' are considered required
        # This is a simplification - adjust based on actual Loculus logic
        required = []
        for field in field_info:
            if field.get("required", True):  # Default to required if not specified
                name = field.get("name")
                if name:
                    required.append(name)
        return required

    def get_file_categories(self) -> List[str]:
        """Get list of file categories supported by this organism.

        Returns:
            List of file category names (e.g., ['nucleotideAlignment', 'siloReads'])
        """
        if self._config is None:
            raise ValueError("Configuration not loaded.")

        schema = self._config.get("schema", {})
        submission_data_types = schema.get("submissionDataTypes", {})
        files_config = submission_data_types.get("files", {})

        if not files_config.get("enabled", False):
            return []

        categories = files_config.get("categories", [])
        if not isinstance(categories, list):
            return []

        return [cat.get("name") for cat in categories if isinstance(cat, dict) and "name" in cat]
