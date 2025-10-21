"""Organism-specific configuration combining schemas and field mappings.

This module provides:
- Unified interface for organism-specific schema management
- Combined validation of SILO schema, Loculus schema, and field mappings
- Easy access to all organism-specific configuration
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

from sr2silo.schema.field_mapping import FieldMapping
from sr2silo.schema.loculus_schema import LoculusSchema
from sr2silo.schema.silo_schema import SiloSchema


class OrganismConfig:
    """Unified organism-specific configuration manager."""

    def __init__(
        self,
        organism: str,
        silo_schema: Optional[SiloSchema] = None,
        loculus_schema: Optional[LoculusSchema] = None,
        field_mapping: Optional[FieldMapping] = None,
    ) -> None:
        """Initialize organism configuration.

        Args:
            organism: Organism identifier (e.g., 'covid', 'sc2')
            silo_schema: SiloSchema instance. If None, must be set later.
            loculus_schema: LoculusSchema instance. If None, must be set later.
            field_mapping: FieldMapping instance. If None, uses default for organism.
        """
        self.organism = organism
        self.silo_schema = silo_schema
        self.loculus_schema = loculus_schema
        self.field_mapping = field_mapping or FieldMapping(organism)

    @classmethod
    def from_api_and_config(
        cls,
        organism: str,
        lapis_url: str,
        loculus_config: Dict[str, Any],
        cache_dir: Optional[Path] = None,
        use_cache: bool = True,
    ) -> "OrganismConfig":
        """Create OrganismConfig by fetching SILO schema from API and loading Loculus config.

        Args:
            organism: Organism identifier
            lapis_url: Base URL for the LAPIS API
            loculus_config: Loculus configuration dictionary (organism-specific portion)
            cache_dir: Optional directory for caching SILO schema
            use_cache: If True, use cached SILO schema when available

        Returns:
            OrganismConfig instance with loaded schemas
        """
        # Initialize SILO schema and fetch from API
        silo_schema = SiloSchema(lapis_url, organism, cache_dir)
        silo_schema.fetch_schema(use_cache=use_cache)

        # Load Loculus schema from config dict
        loculus_schema = LoculusSchema.from_dict(organism, loculus_config)

        # Create field mapping with default for organism
        field_mapping = FieldMapping(organism)

        return cls(
            organism=organism,
            silo_schema=silo_schema,
            loculus_schema=loculus_schema,
            field_mapping=field_mapping,
        )

    @classmethod
    def from_api_and_file(
        cls,
        organism: str,
        lapis_url: str,
        loculus_config_file: Path,
        cache_dir: Optional[Path] = None,
        use_cache: bool = True,
    ) -> "OrganismConfig":
        """Create OrganismConfig by fetching SILO schema from API and loading from file.

        Args:
            organism: Organism identifier
            lapis_url: Base URL for the LAPIS API
            loculus_config_file: Path to Loculus YAML configuration file
            cache_dir: Optional directory for caching SILO schema
            use_cache: If True, use cached SILO schema when available

        Returns:
            OrganismConfig instance with loaded schemas
        """
        # Initialize SILO schema and fetch from API
        silo_schema = SiloSchema(lapis_url, organism, cache_dir)
        silo_schema.fetch_schema(use_cache=use_cache)

        # Load Loculus schema from file
        loculus_schema = LoculusSchema.from_file(organism, loculus_config_file)

        # Create field mapping with default for organism
        field_mapping = FieldMapping(organism)

        return cls(
            organism=organism,
            silo_schema=silo_schema,
            loculus_schema=loculus_schema,
            field_mapping=field_mapping,
        )

    def validate(self, strict: bool = False) -> Dict[str, Any]:
        """Validate the complete organism configuration.

        Checks:
        - SILO schema is loaded
        - Loculus schema is loaded
        - Field mapping is valid against both schemas
        - All required fields are mapped

        Args:
            strict: If True, performs additional strict validation

        Returns:
            Dictionary with validation results including:
            - 'valid': Overall validation status
            - 'errors': List of error messages
            - 'warnings': List of warning messages
            - 'field_mapping_validation': Detailed field mapping validation results
        """
        errors = []
        warnings = []

        # Check if schemas are loaded
        if self.silo_schema is None or self.silo_schema._schema is None:
            errors.append("SILO schema not loaded")
        if self.loculus_schema is None or self.loculus_schema._config is None:
            errors.append("Loculus schema not loaded")

        # If schemas are not loaded, cannot continue validation
        if errors:
            return {
                "valid": False,
                "errors": errors,
                "warnings": warnings,
                "field_mapping_validation": {},
            }

        # Get field lists from schemas
        silo_metadata_fields = self.silo_schema.get_metadata_fields()
        loculus_metadata_fields = self.loculus_schema.get_metadata_fields()

        # Validate field mapping
        mapping_validation = self.field_mapping.validate_against_schemas(
            silo_metadata_fields, loculus_metadata_fields, strict=strict
        )

        # Check for critical mapping errors
        if mapping_validation["missing_silo"]:
            errors.append(
                f"Mapped fields not in SILO schema: {mapping_validation['missing_silo']}"
            )
        if mapping_validation["missing_loculus"]:
            errors.append(
                f"Mapped fields not in Loculus schema: {mapping_validation['missing_loculus']}"
            )

        # Check if all Loculus required fields are covered by mapping
        # (This ensures we can submit valid data to Loculus)
        loculus_required = set(self.loculus_schema.get_required_fields())
        mapped_loculus = set(self.field_mapping.get_loculus_fields())
        missing_required = loculus_required - mapped_loculus

        if missing_required:
            # Filter out system fields that are auto-generated (not from SILO)
            system_fields = {"date", "submissionId"}  # These are added during submission
            user_missing = missing_required - system_fields
            if user_missing:
                warnings.append(
                    f"Required Loculus fields not in mapping: {list(user_missing)}"
                )

        # Strict validation warnings
        if strict:
            if mapping_validation["unmapped_loculus"]:
                # Loculus fields that aren't mapped - might want to populate them
                warnings.append(
                    f"Loculus fields not mapped from SILO: {mapping_validation['unmapped_loculus']}"
                )

        # Overall validation status
        is_valid = len(errors) == 0

        return {
            "valid": is_valid,
            "errors": errors,
            "warnings": warnings,
            "field_mapping_validation": mapping_validation,
        }

    def get_metadata_for_submission(self, silo_metadata: Dict[str, Any]) -> Dict[str, Any]:
        """Transform SILO metadata to Loculus format using field mapping.

        Args:
            silo_metadata: Dictionary of metadata from SILO (with snake_case keys)

        Returns:
            Dictionary of metadata for Loculus submission (with camelCase keys)
        """
        loculus_metadata = {}

        for silo_field, loculus_field in self.field_mapping._mapping.items():
            if silo_field in silo_metadata:
                value = silo_metadata[silo_field]
                loculus_metadata[loculus_field] = value
            else:
                logging.warning(
                    f"SILO field '{silo_field}' in mapping but not in provided metadata"
                )

        return loculus_metadata

    def log_configuration_summary(self) -> None:
        """Log a summary of the organism configuration."""
        logging.info(f"=== Organism Configuration Summary: {self.organism} ===")

        if self.silo_schema and self.silo_schema._schema:
            silo_fields = self.silo_schema.get_metadata_fields()
            logging.info(f"SILO metadata fields ({len(silo_fields)}): {silo_fields}")
        else:
            logging.warning("SILO schema not loaded")

        if self.loculus_schema and self.loculus_schema._config:
            loculus_fields = self.loculus_schema.get_metadata_fields()
            logging.info(f"Loculus metadata fields ({len(loculus_fields)}): {loculus_fields}")
        else:
            logging.warning("Loculus schema not loaded")

        mapped_silo = self.field_mapping.get_silo_fields()
        mapped_loculus = self.field_mapping.get_loculus_fields()
        logging.info(f"Field mappings ({len(mapped_silo)}): {dict(zip(mapped_silo, mapped_loculus))}")

        # Run validation and log results
        validation = self.validate(strict=False)
        if validation["valid"]:
            logging.info("✓ Configuration is valid")
        else:
            logging.error(f"✗ Configuration has errors: {validation['errors']}")

        if validation["warnings"]:
            logging.warning(f"Warnings: {validation['warnings']}")
