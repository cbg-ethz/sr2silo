"""Organism-specific field mapping between SILO and Loculus.

This module handles:
- Defining explicit mappings from SILO fields to Loculus fields
- Validating field mappings for completeness
- Ensuring all mapped fields exist in both schemas
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional


class FieldMapping:
    """Manages field mapping between SILO and Loculus for a specific organism."""

    def __init__(self, organism: str, mapping: Optional[Dict[str, str]] = None) -> None:
        """Initialize field mapping.

        Args:
            organism: Organism identifier (e.g., 'covid', 'sc2')
            mapping: Dictionary mapping SILO field names to Loculus field names.
                     If None, uses default mapping for the organism.
        """
        self.organism = organism
        self._mapping = mapping or self._get_default_mapping(organism)

    @staticmethod
    def _get_default_mapping(organism: str) -> Dict[str, str]:
        """Get default field mapping for an organism.

        Args:
            organism: Organism identifier

        Returns:
            Dictionary mapping SILO field names (snake_case) to Loculus field names (camelCase)
        """
        # Default mapping for COVID/SARS-CoV-2
        # Based on the existing code in submit_to_loculus.py
        if organism in ["covid", "sc2", "sars-cov-2"]:
            return {
                # SILO field (snake_case) -> Loculus field (camelCase)
                "sample_id": "sampleId",
                "batch_id": "batchId",
                "location_code": "locationCode",
                "location_name": "locationName",
                "sampling_date": "samplingDate",
                "sr2silo_version": "sr2siloVersion",
            }
        else:
            logging.warning(
                f"No default mapping defined for organism '{organism}'. "
                "Using empty mapping. Please define explicit mapping."
            )
            return {}

    def silo_to_loculus(self, silo_field: str) -> Optional[str]:
        """Convert SILO field name to Loculus field name.

        Args:
            silo_field: SILO field name (typically snake_case)

        Returns:
            Loculus field name (typically camelCase) or None if not mapped
        """
        return self._mapping.get(silo_field)

    def loculus_to_silo(self, loculus_field: str) -> Optional[str]:
        """Convert Loculus field name to SILO field name.

        Args:
            loculus_field: Loculus field name (typically camelCase)

        Returns:
            SILO field name (typically snake_case) or None if not mapped
        """
        reverse_mapping = {v: k for k, v in self._mapping.items()}
        return reverse_mapping.get(loculus_field)

    def get_silo_fields(self) -> List[str]:
        """Get list of all SILO fields in the mapping.

        Returns:
            List of SILO field names
        """
        return list(self._mapping.keys())

    def get_loculus_fields(self) -> List[str]:
        """Get list of all Loculus fields in the mapping.

        Returns:
            List of Loculus field names
        """
        return list(self._mapping.values())

    def add_mapping(self, silo_field: str, loculus_field: str) -> None:
        """Add or update a field mapping.

        Args:
            silo_field: SILO field name
            loculus_field: Loculus field name
        """
        self._mapping[silo_field] = loculus_field
        logging.debug(f"Added mapping: {silo_field} -> {loculus_field}")

    def remove_mapping(self, silo_field: str) -> None:
        """Remove a field mapping.

        Args:
            silo_field: SILO field name to remove
        """
        if silo_field in self._mapping:
            del self._mapping[silo_field]
            logging.debug(f"Removed mapping for: {silo_field}")

    def validate_against_schemas(
        self,
        silo_fields: List[str],
        loculus_fields: List[str],
        strict: bool = True,
    ) -> Dict[str, List[str]]:
        """Validate field mapping against actual schema fields.

        Args:
            silo_fields: List of field names available in SILO schema
            loculus_fields: List of field names available in Loculus schema
            strict: If True, reports warnings for unmapped fields in either schema

        Returns:
            Dictionary with validation results:
            - 'missing_silo': SILO fields in mapping but not in schema
            - 'missing_loculus': Loculus fields in mapping but not in schema
            - 'unmapped_silo': SILO fields in schema but not in mapping (if strict)
            - 'unmapped_loculus': Loculus fields in schema but not in mapping (if strict)
        """
        silo_set = set(silo_fields)
        loculus_set = set(loculus_fields)
        mapped_silo = set(self.get_silo_fields())
        mapped_loculus = set(self.get_loculus_fields())

        results = {
            "missing_silo": list(mapped_silo - silo_set),
            "missing_loculus": list(mapped_loculus - loculus_set),
            "unmapped_silo": [],
            "unmapped_loculus": [],
        }

        if strict:
            # Find fields in schema but not in mapping
            # Note: Not all SILO fields need to be mapped (e.g., read_id, sequences)
            # Only metadata fields that should be submitted to Loculus need mapping
            results["unmapped_silo"] = list(silo_set - mapped_silo)
            results["unmapped_loculus"] = list(loculus_set - mapped_loculus)

        # Log validation results
        if results["missing_silo"]:
            logging.error(
                f"Fields in mapping not found in SILO schema: {results['missing_silo']}"
            )
        if results["missing_loculus"]:
            logging.error(
                f"Fields in mapping not found in Loculus schema: {results['missing_loculus']}"
            )
        if strict:
            if results["unmapped_silo"]:
                logging.info(
                    f"SILO fields not in mapping: {results['unmapped_silo']}"
                )
            if results["unmapped_loculus"]:
                logging.warning(
                    f"Loculus fields not in mapping: {results['unmapped_loculus']}"
                )

        return results

    def is_valid(self, silo_fields: List[str], loculus_fields: List[str]) -> bool:
        """Check if the mapping is valid against schemas.

        Args:
            silo_fields: List of field names available in SILO schema
            loculus_fields: List of field names available in Loculus schema

        Returns:
            True if all mapped fields exist in their respective schemas
        """
        validation = self.validate_against_schemas(
            silo_fields, loculus_fields, strict=False
        )
        return not validation["missing_silo"] and not validation["missing_loculus"]
