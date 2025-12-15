"""Converts V-PIPE's outputs of nucleotide read sequences output to
ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import json
import logging
import os
import tempfile
from pathlib import Path

import pysam

from sr2silo.process import (
    bam_to_sam,
    paired_end_read_merger,
    parse_translate_align_in_batches,
    sam_to_bam,
    sort_bam_file,
)
from sr2silo.vpipe import Sample


def count_bam_reads(bam_file: Path) -> int:
    """Count the number of reads in a BAM file.

    Uses pysam.AlignmentFile.count() which is equivalent to `samtools view -c`.

    Args:
        bam_file: Path to the BAM file.

    Returns:
        Number of reads in the BAM file.
    """
    try:
        with pysam.AlignmentFile(str(bam_file), "rb") as bam:
            return bam.count(until_eof=True)
    except Exception as e:
        logging.warning(f"Could not count reads in {bam_file}: {e}")
        return 0


def nuc_align_to_silo_njson(
    input_file: Path,
    sample_id: str,
    timeline_file: Path,
    output_fp: Path,
    nuc_ref_fp: Path,
    aa_ref_fp: Path,
    skip_merge: bool = False,
    version_info: str | None = None,
    organism: str = "covid",
) -> bool:
    """Process a given input file.

    Args:
        input_file (Path): The file to process.
        sample_id (str): Sample ID to use for metadata.
        timeline_file (Path): The timeline file to cross-reference the metadata.
        output_fp (Path): Path to the output file.
        nuc_ref_fp (str): Filepath to the nucleotide reference i.e. .fasta.
        aa_ref_fp (str): Filepath to the amino acid reference i.e. .fasta.
        skip_merge (bool): Whether to skip merging of paired-end reads.
                           Default is False.
        version_info (str | None): Version information to include in metadata.
                           Default is None.
        organism (str): The organism identifier (e.g., 'covid', 'rsva').
                       Used for timeline column mappings. Default is 'covid'.

    Returns:
        bool: True if processing was performed, False if skipped (0 reads).
              Writes results to the output file (empty if skipped).
    """
    logging.info(f"Current working directory: {os.getcwd()}")

    # check that the file exists
    if not input_file.exists():
        logging.error(f"Input file not found: {input_file}")
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Early check: count reads in input BAM
    read_count = count_bam_reads(input_file)
    logging.info(f"Input BAM contains {read_count} aligned reads")

    if read_count == 0:
        logging.warning(
            f"⏭️  SKIPPED: Sample {sample_id} has 0 aligned reads in input BAM"
        )
        logging.warning(f"   Input file: {input_file}")
        logging.warning("   This is expected for samples with low/no viral load.")
        logging.warning("   Creating empty output file.")
        # Create empty output file so Snakemake workflow completes
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        output_fp.touch()
        return False

    # get the result directory
    result_dir = input_file.parent / "results"
    result_dir.mkdir(parents=True, exist_ok=True)
    # check that output_fp ends with .ndjson.zst
    if output_fp.suffixes != [".ndjson", ".zst"]:
        logging.warning(
            f"Output file extension changed from {output_fp.suffix} to .ndjson.zst"
        )
        output_fp = output_fp.with_suffix(".ndjson.zst")

    logging.info(f"Processing file: {input_file}")

    ##### Get Sample and Batch metadata and write to a file #####
    sample_to_process = Sample(sample_id, organism=organism)
    sample_to_process.enrich_metadata(timeline_file)
    metadata = sample_to_process.get_metadata()

    # Add version information to metadata if provided
    if version_info is not None:
        metadata["sr2silo_version"] = version_info
        logging.info(f"Added version information to metadata: {version_info}")

    metadata_file = result_dir / "metadata.json"
    result_dir.mkdir(parents=True, exist_ok=True)
    with metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)
    logging.info(f"Metadata saved to: {metadata_file}")

    if skip_merge:
        logging.info("Read pair merging step skipped as requested.")
        # Use the input file directly,
        # but make sure it's converted to BAM format if needed
        if input_file.suffix.lower() == ".sam":
            merged_reads_fp = result_dir / f"{input_file.stem}.bam"
            logging.debug("Converting input SAM to BAM format")
            sam_to_bam(input_file, merged_reads_fp)
        else:
            # If it's already a BAM file, just use it
            merged_reads_fp = input_file
    else:
        #####  Merge & Pair reads #####
        merged_reads_sam_fp = result_dir / f"{input_file.stem}_merged.sam"

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir = Path(tmp_dir)
            logging.info("=== Merging and pairing reads using temporary directory ===")

            logging.debug("Sort by QNAME for matching")
            input_bam_sorted_by_qname_fp = (
                tmp_dir / f"{input_file.stem}.sorted_by_qname{input_file.suffix}"
            )
            sort_bam_file(input_file, input_bam_sorted_by_qname_fp, sort_by_qname=True)

            logging.debug("Decompressing input file to SAM")
            input_sam_fp = tmp_dir / f"{input_file.stem}.sam"
            bam_to_sam(input_bam_sorted_by_qname_fp, input_sam_fp)
            logging.debug(f"Decompressed reads saved to: {input_sam_fp}")

            logging.debug("Starting to merge paired-end reads")
            paired_end_read_merger(
                nuc_align_sam_fp=input_sam_fp,
                ref_genome_fasta_fp=nuc_ref_fp,
                output_merged_sam_fp=merged_reads_sam_fp,
            )
            logging.debug(
                f"Merged reads saved to temporary file: {merged_reads_sam_fp}"
            )

        logging.debug("Re-Compressing merged reads to BAM")
        merged_reads_fp = merged_reads_sam_fp.with_suffix(".bam")
        sam_to_bam(merged_reads_sam_fp, merged_reads_fp)
        merged_reads_sam_fp.unlink()
        logging.info(f"Re-Compressed reads saved to: {merged_reads_fp}")

    ##### Translate / Align / Normalize to JSON #####
    logging.info("=== Start translating, aligning and normalizing reads to JSON ===")
    aligned_reads_fp = output_fp
    try:
        aligned_reads_fp = parse_translate_align_in_batches(
            nuc_reference_fp=nuc_ref_fp,
            aa_reference_fp=aa_ref_fp,
            nuc_alignment_fp=merged_reads_fp,
            metadata_fp=metadata_file,
            output_fp=aligned_reads_fp,
        )
    finally:
        # Only remove temporary files if we created them (i.e., if we merged the reads)
        if (
            not skip_merge
            and merged_reads_fp.exists()
            and merged_reads_fp != input_file
        ):
            merged_reads_fp.unlink()
            logging.info(f"Temporary file {merged_reads_fp} removed.")

    logging.info(f"Processed reads saved to: {aligned_reads_fp}")
    logging.info(
        "Processing completed. Use 'submit-to-loculus' command to upload "
        "and submit to SILO."
    )
    return True
