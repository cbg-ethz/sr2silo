"""Entry point for the sr2silo CLI."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Annotated

import typer

from sr2silo.config import (
    get_primer_file,
    get_reference,
    get_timeline_file,
    get_version,
    is_ci_environment,
)
from sr2silo.process_from_vpipe import nuc_align_to_silo_njson
from sr2silo.submit_to_loculus import submit_to_silo, upload_to_s3

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
def run():
    """
    Wrangel short-reads into cleartext alignments,
    optionally translation and align in amino acids.
    """
    typer.echo("Not yet implemented.")


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
    batch_id: Annotated[
        str,
        typer.Option(
            "--batch-id",
            "-b",
            help="Batch ID to use for metadata.",
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
        Path | None,
        typer.Option(
            "--timeline-file",
            "-t",
            help=(
                "Path to the timeline file. "
                "Can be set via TIMELINE_FILE environment variable."
            ),
        ),
    ] = None,
    primer_file: Annotated[
        Path | None,
        typer.Option(
            "--primer-file",
            "-p",
            help=(
                "Path to the primers file. "
                "Can be set via PRIMER_FILE environment variable."
            ),
        ),
    ] = None,
    reference: Annotated[
        str | None,
        typer.Option(
            "--reference",
            "-r",
            help=(
                "See folder names in resources/. "
                "Can be set via NEXTCLADE_REFERENCE environment variable."
            ),
        ),
    ] = None,
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

    # Handle environment variable defaults
    if timeline_file is None:
        env_timeline = get_timeline_file()
        if env_timeline is not None:
            timeline_file = Path(env_timeline)
        else:
            typer.echo(
                "Error: --timeline-file is required or set TIMELINE_FILE "
                "environment variable"
            )
            raise typer.Exit(1)

    if primer_file is None:
        env_primer = get_primer_file()
        if env_primer is not None:
            primer_file = Path(env_primer)
        else:
            typer.echo(
                "Error: --primer-file is required or set PRIMER_FILE "
                "environment variable"
            )
            raise typer.Exit(1)

    if reference is None:
        reference = get_reference()  # This has a default of "sars-cov-2"

    logging.info(f"Processing input file: {input_file}")
    logging.info(f"Using timeline file: {timeline_file}")
    logging.info(f"Using primers file: {primer_file}")
    logging.info(f"Using output file: {output_fp}")
    logging.info(f"Using genome reference: {reference}")
    logging.info(f"Using sample_id: {sample_id}")
    logging.info(f"Using batch_id: {batch_id}")
    logging.info(f"Skip read pair merging: {skip_merge}")

    # check if $TMPDIR is set, if not use /tmp
    if "TMPDIR" in os.environ:
        temp_dir = Path(os.environ["TMPDIR"])
        logging.info(f"Recognize temporary directory set in Env: {temp_dir}")
        logging.info(
            "This will be used for amino acid translation and alignment - by diamond."
        )

    ci_env = is_ci_environment()
    logging.info(f"Running in CI environment: {ci_env}")

    # Get version information if needed
    version_info = get_version(True)

    logging.info(f"Running version: {version_info}")

    nuc_align_to_silo_njson(
        input_file=input_file,
        sample_id=sample_id,
        batch_id=batch_id,
        timeline_file=timeline_file,
        primers_file=primer_file,
        output_fp=output_fp,
        reference=reference,
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
    sample_id: Annotated[
        str,
        typer.Option(
            "--sample-id",
            "-s",
            help="Sample ID for the processed file.",
        ),
    ],
) -> None:
    """
    Upload processed file to S3 and submit to SILO/Loculus.
    """
    typer.echo("Starting upload and submission to SILO.")

    logging.info(f"Processing file: {processed_file}")
    logging.info(f"Using sample_id: {sample_id}")

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

    # Upload to S3 and submit to SILO
    s3_link = upload_to_s3(processed_file, sample_id)
    success = submit_to_silo(result_dir, s3_link)

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
