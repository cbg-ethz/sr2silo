"""Converts V-PIPE's outputs of nucleotide read sequencs output to
   ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import click
import yaml

from sr2silo.config import is_ci_environment
from sr2silo.process import parse_translate_align
from sr2silo.s3 import compress_bz2, upload_file_to_s3
from sr2silo.silo import LapisClient, Submission
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


def process_directory(
    input_dir: Path,
    sample_id: str,
    batch_id: str,
    result_dir: Path,
    nuc_reference: str,
    aa_reference: str,
    timeline_file: Path,
    primers_file: Path,
    file_name: str = "REF_aln_trim.bam",
    database_config: Path = Path("scripts/database_config.yaml"),
) -> None:
    """Process all files in a given directory.

    Args:
        input_dir (Path): The directory to process. i.e. the directory containing the BAM file.
                          to reach samples/A1_05_2024_10_08/20241024_2411515907/alignments/
        result_dir (Path): The directory to save the results.
        nuc_reference (str): The nucliotide reference from the resources folder.
                          see resources/ e.g. "sars-cov-2"
        aa_reference (str): The amino acid reference from the resources folder.
                            see resources/ e.g. "sars-cov-2"
        timeline_file (Path): The timeline file to cross-reference the metadata.
        primers_file (Path): The primers file to cross-reference the metadata.
        file_name (str): The name of the file to process

    Returns:
        None (writes results to the result_dir)
    """

    # TODO: absolb all these intermediary files into a temporary directory

    ## print PWD
    logging.info(f"Current working directory: {os.getcwd()}")

    # check that one was given a directory and not a file and it exists
    if not input_dir.is_dir():
        logging.error(f"Input directory not found, is it a directory?: {input_dir}")
        raise FileNotFoundError(f"Directory not found: {input_dir}")

    logging.info(f"Processing directory: {input_dir}")
    logging.info(f"Assuming the input file is: {file_name}")
    # check that the file exists and also it's .bai file
    sample_fp = input_dir / file_name
    if not sample_fp.exists():
        logging.error(f"Input file not found: {sample_fp}")
        raise FileNotFoundError(f"Input file not found: {sample_fp}")

    ##### Get Sample and Batch metadata and write to a file #####
    sample_to_process = Sample(sample_id, batch_id)
    sample_to_process.enrich_metadata(timeline_file, primers_file)
    metadata = sample_to_process.get_metadata()
    # add nextclade reference to metadata
    resource_fp = Path("./resources") / nuc_reference
    nuc_reference_fp = resource_fp / "nuc_reference_genomes.fasta"
    aa_reference_fp = resource_fp / "aa_reference_genomes.fasta"

    metadata["nuc_reference"] = nuc_reference
    metadata["aa_reference"] = aa_reference
    metadata_file = result_dir / "metadata.json"
    result_dir.mkdir(parents=True, exist_ok=True)
    with metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)
    logging.info(f"Metadata saved to: {metadata_file}")

    #####  Merge & Pair reads #####

    ## TODO: to implement from smallgenomeutils

    ##### Translate / Align / Normalize to JSON #####

    aligned_reads = parse_translate_align(nuc_reference_fp, aa_reference_fp, sample_fp)

    # TODO wrangle the aligned reads to aligned_reads_with_metadata and write to a file

    # write the aligned reads to a file
    aligned_reads_fp = result_dir / "silo_input.ndjson"

    with aligned_reads_fp.open("w") as f:
        for read in aligned_reads.values():
            try:
                f.write(read.to_silo_json() + "\n")
            except Exception as e:
                logging.error(f"Error writing read to file SILO JSON {e}")
                logging.error(f"Read ID: {read.read_id}")
                logging.error(f"Read: {read}")

    #####   Compress & Upload to S3  #####
    file_to_upload = aligned_reads_fp
    compressed_file = result_dir / "silo_input.ndjson.bz2"
    logging.info(f"Compressing file: {file_to_upload}")
    compress_bz2(file_to_upload, compressed_file)

    #  Upload as generate a file name for the submission file, i.e. use the SAMPLE_ID
    logging.info(f"Uploading to S3: {compressed_file}")
    s3_file_name = f"{sample_id}.ndjson.bz2"
    s3_bucket = "sr2silo01"
    s3_link = f"s3://{s3_bucket}/{s3_file_name}"
    upload_file_to_s3(compressed_file, s3_bucket, s3_file_name)

    ##### Submit S3 reference to SILO #####
    logging.info(f"Submitting to Loculus")
    input_fp = make_submission_file(result_dir, s3_link)

    if is_ci_environment():
        logging.info(
            "Running in CI environment, mocking S3 upload, skipping LAPIS submission."
        )
        # set some mock environment variables for Keycloak and submission URLs
        KEYCLOAK_TOKEN_URL = "https://authentication-wise-seqs.loculus.org/realms/loculus/protocol/openid-connect/token"
        SUBMISSION_URL = "https://backend-wise-seqs.loculus.org/test/submit?groupId={group_id}&dataUseTermsType=OPEN"
    else:
        if os.getenv("KEYCLOAK_TOKEN_URL") or os.getenv("SUBMISSION_URL"):
            KEYCLOAK_TOKEN_URL = os.getenv("KEYCLOAK_TOKEN_URL")
            SUBMISSION_URL = os.getenv("SUBMISSION_URL")
        else:
            logging.warning("KEYCLOAK_TOKEN_URL and SUBMISSION_URL not set.")
            logging.warning("Using default values.")
            KEYCLOAK_TOKEN_URL = "https://authentication-wise-seqs.loculus.org/realms/loculus/protocol/openid-connect/token"
            SUBMISSION_URL = "https://backend-wise-seqs.loculus.org/test/submit?groupId={group_id}&dataUseTermsType=OPEN"

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
    else:
        logging.error(f"Error submitting data to Lapis: {response}")
        logging.error(f"Response: {response.text}")


@click.command()
@click.option("--sample_dir", envvar="SAMPLE_DIR", help="Path to the sample directory.")
@click.option("--sample_id", envvar="SAMPLE_ID", help="sample_id to use for metadata.")
@click.option("--batch_id", envvar="BATCH_ID", help="batch_id to use for metadata.")
@click.option(
    "--result_dir", envvar="RESULTS_DIR", help="Path to the results directory."
)
@click.option(
    "--timeline_file", envvar="TIMELINE_FILE", help="Path to the timeline file."
)
@click.option("--primer_file", envvar="PRIMER_FILE", help="Path to the primers file.")
@click.option(
    "--nuc_reference",
    envvar="NUC_REFERENCE",
    default="sars-cov-2",
    help="see folder names in resources/",
)
@click.option(
    "--aa_reference",
    envvar="AA_REFERENCE",
    default="sars-cov-2",
    help="see folder names in resources/",
)
def main(
    sample_dir,
    sample_id,
    batch_id,
    result_dir,
    timeline_file,
    primer_file,
    nuc_reference,
    aa_reference,
    database_config: Path = Path("scripts/database_config.yaml"),
):
    """Process a sample directory."""
    logging.info(f"Processing sample directory: {sample_dir}")
    logging.info(f"Saving results to: {result_dir}")
    logging.info(f"Using timeline file: {timeline_file}")
    logging.info(f"Using primers file: {primer_file}")
    logging.info(f"Using nucliotide reference: {nuc_reference}")
    logging.info(f"Using amino acid reference: {aa_reference}")
    logging.info(f"Using sample_id: {sample_id}")
    logging.info(f"Using batch_id: {batch_id}")
    logging.info(f"Using database_config: {database_config}")

    ci_env = is_ci_environment()
    logging.info(f"Running in CI environment: {ci_env}")

    process_directory(
        input_dir=Path(sample_dir),
        sample_id=sample_id,
        batch_id=batch_id,
        result_dir=Path(result_dir),
        timeline_file=Path(timeline_file),
        primers_file=Path(primer_file),
        nuc_reference=nuc_reference,
        aa_reference=aa_reference,
        database_config=Path(database_config),
    )


if __name__ == "__main__":
    main()
