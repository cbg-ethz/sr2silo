"""
Fixtures for the process module.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Dict

import pytest

import sr2silo.process.translate_align as translate_align
from sr2silo.process.interface import AlignedRead


@pytest.fixture
def aligned_reads() -> Dict[str, AlignedRead]:
    """
    Small mock data with 42 real reads from the combined.bam file.

    Current dataset misses Amino Acid Insertions - i.e. not tested here.
    """

    nuc_ref_fp = Path("resources/sars-cov-2/nuc_reference_genomes.fasta")
    aa_ref_fp = Path("resources/sars-cov-2/aa_reference_genomes.fasta")
    nuc_alignment_fp = Path("tests/data/bam/combined.bam")

    aligned_reads = translate_align.parse_translate_align(
        nuc_ref_fp, aa_ref_fp, nuc_alignment_fp
    )

    from sr2silo.process import enrich_read_with_metadata

    # make mock metadata with all empty strings, but readId
    # print emoptry ReadMetadata scehma to json
    metadata = {
        "read_id": "readId",
        "sample_id": "",
        "batch_id": "",
        "sampling_date": "",
        "sequencing_date": "",
        "location_name": "",
        "read_length": "",
        "primer_protocol": "",
        "location_code": "",
        "flow_cell_serial_number": "",
        "sequencing_well_position": "",
        "primer_protocol_name": "",
        "nextclade_reference": "",
    }

    metadata_fp = Path("tests/data/process/metadata.json")

    # create a temp file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
        json.dump(metadata, f, indent=4)
        metadata_fp = Path(f.name)

    aligned_reads = enrich_read_with_metadata(aligned_reads, metadata_fp)

    return aligned_reads
