"""Entry point for the sr2silo CLI."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Annotated

import typer

from sr2silo.config import (
    get_group_id,
    get_keycloak_token_url,
    get_password,
    get_submission_url,
    get_timeline_file,
    get_username,
    get_version,
    is_ci_environment,
)
from sr2silo.loculus.lapis import LapisClient
from sr2silo.process_from_vpipe import nuc_align_to_silo_njson
from sr2silo.submit_to_loculus import submit_to_silo

# Use force=True to override any existing logging configuration
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", force=True
)


app = typer.Typer(
    name="sr2silo",
    help=(
        "Convert Short-Read nucleotide .bam alignments to cleartext alignments, "
        "with amino acids and insertions, in JSON format."
    ),
    no_args_is_help=False,  # Changed to False so our callback handles no args
)


@app.callback(invoke_without_command=True)
def callback(ctx: typer.Context):
    """Callback function that runs when no subcommand is provided."""
    if ctx.invoked_subcommand is None:
        typer.echo("Well, you gotta decide what to do.. see --help for subcommands")


@app.command()
def process_from_vpipe(
    input_file: Annotated[
        Path,
        typer.Option(
            "--input-file",
            "-i",
            help="Path to the input file.",
        ),
    ],
    sample_id: Annotated[
        str,
        typer.Option(
            "--sample-id",
            "-s",
            help="Sample ID to use for metadata.",
        ),
    ],
    output_fp: Annotated[
        Path,
        typer.Option(
            "--output-fp",
            "-o",
            help="Path to the output file. Must end with .ndjson.",
        ),
    ],
    timeline_file: Annotated[
        Path,
        typer.Option(
            "--timeline-file",
            "-t",
            help="Path to the timeline file.",
        ),
    ],
    lapis_url: Annotated[
        str,
        typer.Option(
            "--lapis-url",
            "-r",
            help="URL of LAPIS instance, hosting SILO database."
            "Used to fetch the nucleotide / amino acid reference.",
        ),
    ],
    skip_merge: Annotated[
        bool,
        typer.Option(
            "--skip-merge/--no-skip-merge",
            help="Skip merging of paired-end reads.",
        ),
    ] = False,
) -> None:
    """
    V-PIPE to SILO conversion with amino acids and special metadata.
    Processing only - use 'submit-to-loculus' command to upload and submit to SILO.
    """
    typer.echo("Starting V-PIPE to SILO conversion.")

    # Resolve timeline_file with environment fallback
    if timeline_file is None:
        timeline_file = get_timeline_file()
        if timeline_file is None:
            logging.error(
                "Timeline file must be provided via --timeline-file "
                "or TIMELINE_FILE environment variable"
            )
            raise typer.Exit(1)

    logging.info(f"Processing input file: {input_file}")
    logging.info(f"Using timeline file: {timeline_file}")
    logging.info(f"Using output file: {output_fp}")
    logging.info(f"Using Lapis URL: {lapis_url}")
    logging.info(f"Using sample_id: {sample_id}")
    logging.info(f"Skip read pair merging: {skip_merge}")

    # check if $TMPDIR is set, if not use /tmp
    if "TMPDIR" in os.environ:
        temp_dir = Path(os.environ["TMPDIR"])
        logging.info(f"Recognize temporary directory set in Env: {temp_dir}")
        logging.info(
            "This will be used for amino acid translation and "
            "alignment - by diamond."
        )

    ci_env = is_ci_environment()
    logging.info(f"Running in CI environment: {ci_env}")

    # Get version information if needed
    version_info = get_version(True)

    logging.info(f"Running version: {version_info}")

    # if not CI envrionement, get the nucleotide and amino acid references, from Lapis
    if not ci_env:
        # make LapisClient
        lapis = LapisClient(lapis_url)
        # fetch references
        reference = lapis.referenceGenome()
        # convert references to FASTA files
        logging.info("Fetching references from Lapis...")
        # get the domain from the lapis_url
        domain = lapis_url.split("//")[-1].split("/")[0]
        # create the directory if it does not exist
        Path(f"resources/references/{domain}").mkdir(parents=True, exist_ok=True)
        # define the paths for the nucleotide and amino acid references
        nuc_ref_fp = Path(f"resources/references/{domain}/nuc_ref.fasta")
        aa_ref_fp = Path(f"resources/references/{domain}/aa_ref.fasta")
        lapis.referenceGenomeToFasta(
            reference_json_string=json.dumps(reference),
            nucleotide_out_fp=nuc_ref_fp,
            amino_acid_out_fp=aa_ref_fp,
        )
        logging.info(f"Fetched references from Lapis: {nuc_ref_fp} and {aa_ref_fp}")
    else:
        logging.info(
            "Running in CI environment, using default references \
            from resources/references/sars-cov-2/"
        )
        nuc_ref_fp = Path("resources/references/sars-cov-2/nuc_ref.fasta")
        aa_ref_fp = Path("resources/references/sars-cov-2/aa_ref.fasta")

    nuc_align_to_silo_njson(
        input_file=input_file,
        sample_id=sample_id,
        timeline_file=timeline_file,
        output_fp=output_fp,
        nuc_ref_fp=nuc_ref_fp,
        aa_ref_fp=aa_ref_fp,
        skip_merge=skip_merge,
        version_info=version_info,
    )


@app.command()
def submit_to_loculus(
    processed_file: Annotated[
        Path,
        typer.Option(
            "--processed-file",
            "-f",
            help="Path to the processed .ndjson.zst file to upload and submit.",
        ),
    ],
    keycloak_token_url: Annotated[
        str | None,
        typer.Option(
            "--keycloak-token-url",
            help="Keycloak authentication URL. Falls back to "
            "KEYCLOAK_TOKEN_URL environment variable.",
        ),
    ] = None,
    submission_url: Annotated[
        str | None,
        typer.Option(
            "--submission-url",
            help="Loculus submission URL. Falls back to "
            "SUBMISSION_URL environment variable.",
        ),
    ] = None,
    group_id: Annotated[
        int | None,
        typer.Option(
            "--group-id",
            help="Group ID for submission. Falls back to "
            "GROUP_ID environment variable.",
        ),
    ] = None,
    username: Annotated[
        str | None,
        typer.Option(
            "--username",
            help="Username for authentication. Falls back to "
            "USERNAME environment variable.",
        ),
    ] = None,
    password: Annotated[
        str | None,
        typer.Option(
            "--password",
            help="Password for authentication. Falls back to "
            "PASSWORD environment variable.",
        ),
    ] = None,
) -> None:
    """
    Upload processed file to S3 and submit to SILO/Loculus.
    """
    typer.echo("Starting upload and submission to SILO.")

    # Resolve keycloak_token_url with environment fallback
    if keycloak_token_url is None:
        keycloak_token_url = get_keycloak_token_url()

    # Resolve submission_url with environment fallback
    if submission_url is None:
        submission_url = get_submission_url()

    # Resolve group_id with environment fallback
    if group_id is None:
        group_id = get_group_id()

    # Resolve username with environment fallback
    if username is None:
        username = get_username()

    # Resolve password with environment fallback
    if password is None:
        password = get_password()

    logging.info(f"Processing file: {processed_file}")
    logging.info(f"Using Keycloak token URL: {keycloak_token_url}")
    logging.info(f"Using submission URL: {submission_url}")
    logging.info(f"Using group ID: {group_id}")
    logging.info(f"Using username: {username}")

    # Check if the processed file exists
    if not processed_file.exists():
        logging.error(f"Processed file not found: {processed_file}")
        raise typer.Exit(1)

    # Check if file has correct extension
    if processed_file.suffixes != [".ndjson", ".zst"]:
        logging.error(
            f"File must have .ndjson.zst extension, got: {processed_file.suffixes}"
        )
        raise typer.Exit(1)

    ci_env = is_ci_environment()
    logging.info(f"Running in CI environment: {ci_env}")

    # Get version information
    version_info = get_version(True)
    logging.info(f"Running version: {version_info}")

    # Get the result directory (parent of the processed file)
    result_dir = processed_file.parent

    # Submit to SILO using the pre-signed upload approach
    # This will handle both metadata and processed file upload via pre-signed URLs

    success = submit_to_silo(
        result_dir,
        processed_file,
        keycloak_token_url=keycloak_token_url,
        submission_url=submission_url,
        group_id=group_id,
        username=username,
        password=password,
    )

    if success:
        typer.echo("Upload and submission completed successfully.")
    else:
        typer.echo("Upload and submission failed.")
        raise typer.Exit(1)


def main():
    """Main entry point for the sr2silo CLI."""
    typer.echo("Well, you gotta decide what to do.. see --help for subcommands")


if __name__ == "__main__":
    main()
