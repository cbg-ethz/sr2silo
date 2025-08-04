"""SILO-specific pydantic schemas that define the expected format for reads in
    the SILO database.

This module contains the schema definitions used to validate read data before submission
to the SILO database, ensuring all records conform to the expected format.

The schemas follow IUPAC Codes.
https://www.bioinformatics.org/sms2/iupac.html
"""

from __future__ import annotations

import logging
import re
from typing import Dict, List, Optional

from pydantic import BaseModel, RootModel, field_validator, model_validator

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


class ReadMetadata(BaseModel):
    """V-Pipe SILO-specific pydantic schema for ReadMetadata JSON format.
    (specific to WISE / run at ETHZ)
    """

    read_id: str
    sample_id: str
    batch_id: str  # Can be empty string for samples without batch_id
    sampling_date: str
    location_name: str
    read_length: str
    primer_protocol: str
    location_code: str
    sr2silo_version: str


class AlignedNucleotideSequences(BaseModel):
    """SILO-specific pydantic schema for AlignedNucleotideSequences JSON format."""

    main: str


class NucleotideInsertions(BaseModel):
    """SILO-specific pydantic schema for NucleotideInsertions JSON format."""

    main: List[str]

    @field_validator("main")
    @classmethod
    def validate_nuc_insertions_format(cls, v: List[str]) -> List[str]:
        """Validate that nucleotide insertions have the format 'position:sequence'."""
        pattern = r"^\d+:[ACGTN]+$"
        for insertion in v:
            if not re.match(pattern, insertion):
                raise ValueError(
                    f"Nucleotide insertion '{insertion}' is not in the expected "
                    "format. Expected format: 'position : sequence' "
                    "(e.g., '123:ACGT')"
                )
        return v


class AminoAcidSequences(RootModel):
    """SILO-specific pydantic schema for AminoAcidSequences JSON format."""

    root: Dict[str, Optional[str]]


class AminoAcidInsertions(RootModel):
    """SILO-specific pydantic schema for AminoAcidInsertions JSON format."""

    root: Dict[str, List[str]]

    @model_validator(mode="after")
    def validate_aa_insertions_format(self) -> "AminoAcidInsertions":
        """Validate that amino acid insertions have the format 'position:sequence'."""
        # Validate each amino acid insertion using one-letter codes only
        pattern = r"^\d+:[A-Z*\-]+$"
        for gene, insertions in self.root.items():
            for insertion in insertions:
                if not re.match(pattern, insertion):
                    raise ValueError(
                        f"Amino acid insertion '{insertion}' for gene '{gene}' "
                        "is not in the expected format. "
                        "Expected format: 'position:sequence' (e.g., '123:AST')"
                    )
        return self


class NucleotideSegment(BaseModel):
    """Schema for a nucleotide genomic segment (e.g., main nucleotide sequence)."""

    sequence: str
    insertions: List[str]
    offset: int

    @field_validator("sequence")
    @classmethod
    def validate_nucleotide_sequence(cls, v: str) -> str:
        """Validate that sequence contains only valid nucleotide characters."""
        nucleotide_pattern = r"^[ACGTN\-]*$"

        if not re.match(nucleotide_pattern, v, re.IGNORECASE):
            raise ValueError(
                "Nucleotide sequence contains invalid characters. "
                "Expected nucleotides (ACGTN-) only."
            )
        return v

    @field_validator("insertions")
    @classmethod
    def validate_nucleotide_insertions(cls, v: List[str]) -> List[str]:
        """Validate that nucleotide insertions have the format 'position:sequence'."""
        pattern = r"^\d+:[ACGTN]+$"
        for insertion in v:
            if not re.match(pattern, insertion, re.IGNORECASE):
                raise ValueError(
                    f"Nucleotide insertion '{insertion}' is not in the expected "
                    "format. Expected format: 'position:sequence' with nucleotides "
                    "only (e.g., '123:ACGT')"
                )
        return v

    @field_validator("offset")
    @classmethod
    def validate_offset(cls, v: int) -> int:
        """Validate that offset is non-negative."""
        if v < 0:
            raise ValueError("Offset must be non-negative")
        return v


class AminoAcidSegment(BaseModel):
    """Schema for an amino acid genomic segment (e.g., gene sequences)."""

    sequence: str
    insertions: List[str]
    offset: int

    @field_validator("sequence")
    @classmethod
    def validate_amino_acid_sequence(cls, v: str) -> str:
        """Validate that sequence contains only valid amino acid characters."""
        amino_acid_pattern = r"^[A-Z*\-]*$"

        if not re.match(amino_acid_pattern, v):
            raise ValueError(
                "Amino acid sequence contains invalid characters. "
                "Expected amino acids (A-Z*-) only."
            )
        return v

    @field_validator("insertions")
    @classmethod
    def validate_amino_acid_insertions(cls, v: List[str]) -> List[str]:
        """Validate that amino acid insertions have the format 'position:sequence'."""
        pattern = r"^\d+:[A-Z*\-]+$"
        for insertion in v:
            if not re.match(pattern, insertion):
                raise ValueError(
                    f"Amino acid insertion '{insertion}' is not in the expected "
                    "format. Expected format: 'position:sequence' with amino acids "
                    "only (e.g., '45:MYK')"
                )
            # Additional check: ensure the sequence part contains only amino acids,
            # not nucleotides
            position, sequence_part = insertion.split(":", 1)
            if re.match(r"^[ACGTN]+$", sequence_part, re.IGNORECASE):
                raise ValueError(
                    f"Amino acid insertion '{insertion}' contains nucleotides. "
                    "Expected amino acids (A-Z*-) only, not nucleotides (ACGTN)."
                )
        return v

    @field_validator("offset")
    @classmethod
    def validate_offset(cls, v: int) -> int:
        """Validate that offset is non-negative."""
        if v < 0:
            raise ValueError("Offset must be non-negative")
        return v


class AlignedReadSchema(BaseModel):
    """SILO-specific pydantic schema for AlignedRead JSON format.

    This schema validates the new SILO format where:
    - read_id is a required field at the root level
    - Metadata fields are at the root level
    - 'main' is a nucleotide segment with nucleotide-specific validation
    - Gene segments are amino acid segments with amino acid-specific validation
    - Unaligned sequences are prefixed with 'unaligned_'
    """

    # Required read_id field
    read_id: str

    # Main nucleotide segment (required)
    main: NucleotideSegment

    # Additional fields will be validated dynamically
    model_config = {"extra": "allow"}

    @model_validator(mode="after")
    def validate_dynamic_fields(self) -> "AlignedReadSchema":
        """Validate gene segments, unaligned sequences, and metadata fields."""
        # Get ReadMetadata field names for validation
        metadata_fields = set(ReadMetadata.model_fields.keys())

        for field_name, field_value in self.__dict__.items():
            if field_name in ["read_id", "main"]:
                continue

            # Check for gene segments (should be amino acid segments or None)
            if isinstance(field_value, dict) and all(
                k in field_value for k in ["sequence", "insertions", "offset"]
            ):
                # Validate as AminoAcidSegment for gene segments
                try:
                    AminoAcidSegment(**field_value)
                except Exception as e:
                    raise ValueError(f"Invalid amino acid segment '{field_name}': {e}")

            # Allow None for empty gene segments
            elif field_value is None:
                continue

            # Check for unaligned sequences (should be strings with nucleotides)
            elif field_name.startswith("unaligned_") and isinstance(field_value, str):
                # Validate unaligned sequence format (nucleotides only)
                nucleotide_pattern = r"^[ACGTN\-]*$"
                if not re.match(nucleotide_pattern, field_value, re.IGNORECASE):
                    raise ValueError(
                        f"Unaligned sequence '{field_name}' contains invalid "
                        f"characters. Expected nucleotides (ACGTN-) only."
                    )

            # Validate metadata fields if they match ReadMetadata field names
            elif field_name in metadata_fields:
                # Validate individual metadata field types and constraints
                metadata_field_info = ReadMetadata.model_fields[field_name]
                expected_type = metadata_field_info.annotation

                # Basic type checking for string fields
                if expected_type is str and not isinstance(field_value, str):
                    raise ValueError(
                        f"Metadata field '{field_name}' should be a string, "
                        f"got {type(field_value).__name__}"
                    )

            # Allow other metadata fields (strings, numbers)
            elif isinstance(field_value, (str, int, float)):
                continue

            else:
                # For gene segments that aren't properly structured
                if not field_name.startswith("unaligned_") and not isinstance(
                    field_value, (str, int, float)
                ):
                    raise ValueError(
                        f"Field '{field_name}' should be either an amino acid "
                        f"segment (with sequence, insertions, offset), None, or "
                        f"metadata (string/number)"
                    )

        return self
