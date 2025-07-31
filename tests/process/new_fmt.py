from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Dict

import sr2silo.process.translate_align as translate_align
from sr2silo.process.interface import AlignedRead


def aligned_reads(aa_ref_sarscov2_fp, nuc_ref_sarscov2_fp) -> Dict[str, AlignedRead]:
    """
    Small mock data with 42 real reads from the combined.bam file.

    Current dataset misses Amino Acid Insertions - i.e. not tested here.
    """

    nuc_ref_fp = nuc_ref_sarscov2_fp
    aa_ref_fp = aa_ref_sarscov2_fp
    nuc_alignment_fp = Path("../data/bam/combined.bam")

    aligned_reads = translate_align.parse_translate_align(
        nuc_ref_fp, aa_ref_fp, nuc_alignment_fp
    )

    # make mock metadata with all empty strings, but readId
    # print emoptry ReadMetadata scehma to json
    metadata = {
        "read_id": "readId",
        "sample_id": "34",
        "batch_id": "",
        "sampling_date": "",
        "location_name": "",
        "read_length": "",
        "primer_protocol": "",
        "location_code": "",
        "sr2silo_version": "",
    }

    metadata_fp = Path("../data/process/metadata.json")

    # create a temp file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
        json.dump(metadata, f, indent=4)
        metadata_fp = Path(f.name)

    enrich_read_with_metadata = translate_align.curry_read_with_metadata(metadata_fp)

    aligned_reads = {
        read_id: enrich_read_with_metadata(read)
        for read_id, read in aligned_reads.items()
    }

    # Update expected data write to ndjson
    ndjson_fp = Path("../data/process/aligned_reads_flat_meta.ndjson")
    with open(ndjson_fp, "w") as f:
        for read_id, read in aligned_reads.items():
            f.write(f"{read.to_silo_json()}\n")

    return aligned_reads


if __name__ == "__main__":
    # Example usage
    aa_ref_sarscov2_fp = Path("../../resources/references/sars-cov-2/aa_ref.fasta")
    nuc_ref_sarscov2_fp = Path("../../resources/references/sars-cov-2/nuc_ref.fasta")
    reads = aligned_reads(aa_ref_sarscov2_fp, nuc_ref_sarscov2_fp)
    for read_id, read in reads.items():
        print(f"{read_id}")
        print(f"{json.dumps(json.loads(read.to_silo_json()), indent=4)}")
