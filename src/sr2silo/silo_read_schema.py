"""SILO-specific pydantic schemas that define the expected format for reads in
    the SILO database.

This module contains the schema definitions used to validate read data before submission
to the SILO database, ensuring all records conform to the expected format.
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
    batch_id: str
    sampling_date: str
    sequencing_date: str
    location_name: str
    read_length: str
    primer_protocol: str
    location_code: str
    flow_cell_serial_number: str
    sequencing_well_position: str
    primer_protocol_name: str
    nextclade_reference: str


class AlignedNucleotideSequences(BaseModel):
    """SILO-specific pydantic schema for AlignedNucleotideSequences JSON format."""

    main: str


class UnalignedNucleotideSequences(BaseModel):
    """SILO-specific pydantic schema for UnalignedNucleotideSequences JSON format."""

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
        pattern = r"^\d+:[ARNDCEQGHILKMFPSTWYV]+$"
        for gene, insertions in self.root.items():
            for insertion in insertions:
                if not re.match(pattern, insertion):
                    raise ValueError(
                        f"Amino acid insertion '{insertion}' for gene '{gene}'"
                        "is not in the expected format."
                        "Expected format: 'position:sequence' (e.g., '123:AST')"
                    )
        return self


class AlignedReadSchema(BaseModel):
    """SILO-specific pydantic schema for AlignedRead JSON format."""

    metadata: Optional[ReadMetadata] = None
    nucleotideInsertions: NucleotideInsertions
    aminoAcidInsertions: AminoAcidInsertions
    alignedNucleotideSequences: AlignedNucleotideSequences
    unalignedNucleotideSequences: UnalignedNucleotideSequences
    alignedAminoAcidSequences: AminoAcidSequences
