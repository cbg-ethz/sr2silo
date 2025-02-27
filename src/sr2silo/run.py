"""Converts V-PIPE's outputs of nucleotide read sequencs output to
   ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import click

from sr2silo.config import (
    get_keycloak_token_url,
    get_mock_urls,
    get_submission_url,
    is_ci_environment,
)
from sr2silo.process import parse_translate_align_in_batches
from sr2silo.silo import LapisClient, Submission
from sr2silo.storage import upload_file_to_s3
from sr2silo.vpipe import Sample

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def load_config(config_file: Path) -> dict:
    """Load a JSON configuration file."""
    try:
        with config_file.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"Config file not found: {config_file}")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from config file: {config_file} - {e}")
        raise


def make_submission_file(result_dir: Path, srLink: str) -> Path:
    """Create a submission file with the given S3 link.

    Args:
        result_dir (Path): The directory to save the submission file.
        srLink (str): The S3 link to include in the submission file.

    Returns:
        Path: The path to the created submission file.
    """
    result_dir_submission = result_dir / "submission"
    result_dir_submission.mkdir(parents=True, exist_ok=True)

    submission_metadata_fp = result_dir_submission / "metadata.tsv"
    with submission_metadata_fp.open("w") as f:
        f.write("submissionId\ts3Link\tversionComment\n")
        f.write(f"001\t{srLink}\t\n")
    logging.info(f"Submission metadata saved to: {submission_metadata_fp}")

    return submission_metadata_fp


def upload_to_s3(aligned_reads_fp: Path, sample_id: str) -> str:
    """Upload a file to S3 bucket.

    Args:
        aligned_reads_fp (Path): The file to upload.
        sample_id (str): Sample ID used in the S3 filename.

    Returns:
        str: The S3 link to the uploaded file.
    """
    logging.info(f"Uploading to S3: {aligned_reads_fp}")
    suffix = aligned_reads_fp.suffix
    s3_file_name = f"{sample_id}.{suffix}"
    s3_bucket = "sr2silo01"
    s3_link = f"s3://{s3_bucket}/{s3_file_name}"
    upload_file_to_s3(aligned_reads_fp, s3_bucket, s3_file_name)
    return s3_link


def submit_to_silo(result_dir: Path, s3_link: str) -> bool:
    """Submit S3 reference to SILO.

    Args:
        result_dir (Path): Directory where to save submission files.
        s3_link (str): The S3 link to submit.

    Returns:
        bool: True if submission was successful, False otherwise.
    """
    logging.info("Submitting to Loculus")
    input_fp = make_submission_file(result_dir, s3_link)

    if is_ci_environment():
        logging.info(
            "Running in CI environment, mocking S3 upload, skipping LAPIS submission."
        )
        # Get mock URLs for CI environment
        KEYCLOAK_TOKEN_URL, SUBMISSION_URL = get_mock_urls()
    else:
        # Get URLs from environment or use defaults
        KEYCLOAK_TOKEN_URL = get_keycloak_token_url()
        SUBMISSION_URL = get_submission_url()
        logging.info(f"Using Keycloak URL: {KEYCLOAK_TOKEN_URL}")
        logging.info(f"Using submission URL: {SUBMISSION_URL}")

    client = LapisClient(KEYCLOAK_TOKEN_URL, SUBMISSION_URL)  # type: ignore
    client.authenticate(username="testuser", password="testuser")
    submission_ids = Submission.get_submission_ids_from_tsv(input_fp)
    fasta_str = Submission.generate_placeholder_fasta(submission_ids)
    submission = Submission(fasta_str, input_fp)
    response = client.submit(group_id=1, data=submission)
    if response.status_code == 200:
        logging.info("Submission successful.")
        logging.info(
            "You can approve the upload for release at:\n\n"
            "https://wise-seqs.loculus.org/salmonella/submission/1/review"
        )
        return True
    else:
        logging.error(f"Error submitting data to Lapis: {response}")
        logging.error(f"Response: {response.text}")
        return False


def process_file(
    input_file: Path,
    sample_id: str,
    batch_id: str,
    reference: str,
    timeline_file: Path,
    primers_file: Path,
    output_fp: Path,
    upload: bool = False,
) -> None:
    """Process a given input file.

    Args:
        input_file (Path): The file to process.
        sample_id (str): Sample ID to use for metadata.
        batch_id (str): Batch ID to use for metadata.
        reference (str): The nucleotide / amino acid reference from
                         the resources folder.
        timeline_file (Path): The timeline file to cross-reference the metadata.
        primers_file (Path): The primers file to cross-reference the metadata.
        output_fp (Path): Path to the output file.
        upload (bool): Whether to upload and submit to SILO. Default is False.

    Returns:
        None (writes results to the result_dir)
    """
    logging.info(f"Current working directory: {os.getcwd()}")

    # check that the file exists
    if not input_file.exists():
        logging.error(f"Input file not found: {input_file}")
        raise FileNotFoundError(f"Input file not found: {input_file}")

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
    sample_to_process = Sample(sample_id, batch_id)
    sample_to_process.enrich_metadata(timeline_file, primers_file)
    metadata = sample_to_process.get_metadata()
    # add nextclade reference to metadata
    resource_fp = Path("./resources") / reference
    nuc_reference_fp = resource_fp / "nuc_reference_genomes.fasta"
    aa_reference_fp = resource_fp / "aa_reference_genomes.fasta"

    # TODO: get reference from nextclade or loculus
    metadata["nextclade_reference"] = reference
    # metadata["nuc_reference"] = nuc_reference
    # metadata["aa_reference"] = aa_reference
    metadata_file = result_dir / "metadata.json"
    result_dir.mkdir(parents=True, exist_ok=True)
    with metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)
    logging.info(f"Metadata saved to: {metadata_file}")

    #####  Merge & Pair reads #####

    ## TODO: to implement from smallgenomeutils

    ##### Translate / Align / Normalize to JSON #####
    logging.info("Start translating, aligning and normalizing reads to JSON")
    aligned_reads_fp = output_fp
    aligned_reads_fp = parse_translate_align_in_batches(
        nuc_reference_fp=nuc_reference_fp,
        aa_reference_fp=aa_reference_fp,
        nuc_alignment_fp=input_file,
        metadata_fp=metadata_file,
        output_fp=aligned_reads_fp,
    )

    if upload:
        s3_link = upload_to_s3(aligned_reads_fp, sample_id)
        submit_to_silo(result_dir, s3_link)
    else:
        logging.info("Skipping upload and submission to S3 and SILO.")


@click.command()
@click.option("--input_file", envvar="INPUT_FILE", help="Path to the input file.")
@click.option("--sample_id", envvar="SAMPLE_ID", help="sample_id to use for metadata.")
@click.option("--batch_id", envvar="BATCH_ID", help="batch_id to use for metadata.")
@click.option(
    "--timeline_file", envvar="TIMELINE_FILE", help="Path to the timeline file."
)
@click.option("--primer_file", envvar="PRIMER_FILE", help="Path to the primers file.")
@click.option(
    "--output_fp",
    envvar="OUTPUT_FP",
    help="Path to the output file. Must end with .ndjson.",
)
@click.option(
    "--reference",
    envvar="REFERENCE",
    default="sars-cov-2",
    help="see folder names in resources/",
)
@click.option(
    "--upload", is_flag=True, default=False, help="Upload and submit to SILO."
)
def main(
    input_file,
    sample_id,
    batch_id,
    timeline_file,
    primer_file,
    output_fp,
    reference,
    upload,
):
    """Process a sample file."""
    logging.info(f"Processing input file: {input_file}")
    logging.info(f"Using timeline file: {timeline_file}")
    logging.info(f"Using primers file: {primer_file}")
    logging.info(f"Using output file: {output_fp}")
    logging.info(f"Using genome reference: {reference}")
    logging.info(f"Using sample_id: {sample_id}")
    logging.info(f"Using batch_id: {batch_id}")
    logging.info(f"Upload to S3 and submit to SILO: {upload}")

    ci_env = is_ci_environment()
    logging.info(f"Running in CI environment: {ci_env}")

    process_file(
        input_file=Path(input_file),
        sample_id=sample_id,
        batch_id=batch_id,
        timeline_file=Path(timeline_file),
        primers_file=Path(primer_file),
        output_fp=Path(output_fp),
        reference=reference,
        upload=upload,
    )


if __name__ == "__main__":
    main()
