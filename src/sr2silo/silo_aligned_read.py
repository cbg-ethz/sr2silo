"""SILO-specific pydantic schemas for AlignedRead JSON format."""

from __future__ import annotations
import json
import logging
from typing import Any, Dict, List, Optional
from pydantic import BaseModel, RootModel, model_validator

logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")

# --- SILO-specific pydantic schemas for AlignedRead JSON format ---
class ReadMetadata(BaseModel):
    read_id: str
    sequencing_date: str
    location_name: str
    batch_id: str
    read_length: str
    primer_protocol: str
    location_code: str
    flow_cell_serial_number: str
    nextclade_reference: str
    sequencing_well_position: str
    sample_id: str
    sampling_date: str
    primer_protocol_name: str

class AlignedNucleotideSequences(BaseModel):
    main: str

class UnalignedNucleotideSequences(BaseModel):
    main: str

class NucleotideInsertions(BaseModel):
    main: List[str]

class AminoAcidSequences(RootModel):
    root: Dict[str, Optional[str]]

class AminoAcidInsertions(RootModel):
    root: Dict[str, List[str]]

class AlignedReadSchema(BaseModel):
    metadata: Optional[ReadMetadata] = None
    nucleotideInsertions: NucleotideInsertions
    aminoAcidInsertions: AminoAcidInsertions
    alignedNucleotideSequences: AlignedNucleotideSequences
    unalignedNucleotideSequences: UnalignedNucleotideSequences
    alignedAminoAcidSequences: AminoAcidSequences

