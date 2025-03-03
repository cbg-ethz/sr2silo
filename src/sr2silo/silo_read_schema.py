"""SILO-specific pydantic schemas that define the expected format for reads in
    the SILO database.

This module contains the schema definitions used to validate read data before submission
to the SILO database, ensuring all records conform to the expected format.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

from pydantic import BaseModel, RootModel

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


# TODO: note that this is V-Pipe Specific
class ReadMetadata(BaseModel):
    """SILO-specific pydantic schema for ReadMetadata JSON format."""

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


class AminoAcidSequences(RootModel):
    """SILO-specific pydantic schema for AminoAcidSequences JSON format."""

    root: Dict[str, Optional[str]]


class AminoAcidInsertions(RootModel):
    """SILO-specific pydantic schema for AminoAcidInsertions JSON format."""

    root: Dict[str, List[str]]


class AlignedReadSchema(BaseModel):
    """SILO-specific pydantic schema for AlignedRead JSON format."""

    metadata: Optional[ReadMetadata] = None
    nucleotideInsertions: NucleotideInsertions
    aminoAcidInsertions: AminoAcidInsertions
    alignedNucleotideSequences: AlignedNucleotideSequences
    unalignedNucleotideSequences: UnalignedNucleotideSequences
    alignedAminoAcidSequences: AminoAcidSequences
