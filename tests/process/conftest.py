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
def aligned_reads(aa_ref_sarscov2_fp, nuc_ref_sarscov2_fp) -> Dict[str, AlignedRead]:
    """
    Small mock data with 42 real reads from the combined.bam file.

    Current dataset misses Amino Acid Insertions - i.e. not tested here.
    """

    nuc_ref_fp = nuc_ref_sarscov2_fp
    aa_ref_fp = aa_ref_sarscov2_fp
    nuc_alignment_fp = Path("tests/data/bam/combined.bam")

    with tempfile.TemporaryDirectory() as temp_dir:
        aa_db_fp = Path(temp_dir) / "diamond_db.dmnd"
        translate_align._make_diamond_db(aa_ref_fp, aa_db_fp)

        aligned_reads = translate_align.parse_translate_align(
            nuc_ref_fp, aa_ref_fp, nuc_alignment_fp, aa_db_fp
        )

    # make mock metadata with all empty strings, but readId
    # print empty ReadMetadata schema to json
    metadata = {
        "read_id": "readId",
        "sample_id": "",
        "batch_id": "",
        "sampling_date": "",
        "location_name": "",
        "read_length": "",
        "location_code": "",
        "sr2silo_version": "",
    }

    metadata_fp = Path("tests/data/process/metadata.json")

    # create a temp file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
        json.dump(metadata, f, indent=4)
        metadata_fp = Path(f.name)

    enrich_read_with_metadata = translate_align.curry_read_with_metadata(metadata_fp)

    aligned_reads = {
        read_id: enrich_read_with_metadata(read)
        for read_id, read in aligned_reads.items()
    }

    return aligned_reads


@pytest.fixture
def micro_bam_fp() -> Path:
    """Path to the micro BAM file."""
    return Path("tests/data/bam/micro/micro.bam")


@pytest.fixture
def micro_bam_sam_fp() -> Path:
    """Path to the micro BAM file."""
    return Path("tests/data/bam/micro/micro.sam")


@pytest.fixture
def micro_reference_fp() -> Path:
    """Path to the micro reference file."""
    return Path("tests/data/bam/micro/micro_ref.fasta")
