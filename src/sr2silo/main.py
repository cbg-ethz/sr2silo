"""Entry point for the sr2silo CLI."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Annotated

import typer

from sr2silo.config import (
    get_batch_id,
    get_default_input_file,
    get_nextclade_reference,
    get_primer_file,
    get_results_dir,
    get_sample_id,
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
        Path | None,
        typer.Option(
            "--input-file",
            "-i",
            help="Path to the input file. Can be set via SAMPLE_DIR environment variable (will look for REF_aln_trim.bam).",
        ),
    ] = None,
    sample_id: Annotated[
        str | None,
        typer.Option(
            "--sample-id",
            "-s",
            help="Sample ID to use for metadata. Can be set via SAMPLE_ID environment variable.",
        ),
    ] = None,
    batch_id: Annotated[
        str | None,
        typer.Option(
            "--batch-id",
            "-b",
            help="Batch ID to use for metadata. Can be set via BATCH_ID environment variable.",
        ),
    ] = None,
    timeline_file: Annotated[
        Path | None,
        typer.Option(
            "--timeline-file",
            "-t",
            help="Path to the timeline file. Can be set via TIMELINE_FILE environment variable.",
        ),
    ] = None,
    primer_file: Annotated[
        Path | None,
        typer.Option(
            "--primer-file",
            "-p",
            help="Path to the primers file. Can be set via PRIMER_FILE environment variable.",
        ),
    ] = None,
    output_fp: Annotated[
        Path | None,
        typer.Option(
            "--output-fp",
            "-o",
            help="Path to the output file. Must end with .ndjson. Can be set via RESULTS_DIR environment variable (will auto-generate filename).",
        ),
    ] = None,
    reference: Annotated[
        str | None,
        typer.Option(
            "--reference",
            "-r",
            help="See folder names in resources/. Can be set via NEXTCLADE_REFERENCE environment variable.",
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
    # Apply environment variable defaults
    if input_file is None:
        input_file = get_default_input_file()
    if sample_id is None:
        sample_id = get_sample_id()
    if batch_id is None:
        batch_id = get_batch_id()
    if timeline_file is None:
        env_timeline = get_timeline_file()
        if env_timeline:
            timeline_file = Path(env_timeline)
    if primer_file is None:
        env_primer = get_primer_file()
        if env_primer:
            primer_file = Path(env_primer)
    if reference is None:
        reference = get_nextclade_reference()
    if output_fp is None:
        env_results = get_results_dir()
        if env_results and sample_id and batch_id:
            output_fp = Path(env_results) / f"{sample_id}_{batch_id}_silo_input.ndjson"

    # Validate required parameters
    missing_params = []
    if input_file is None:
        missing_params.append("--input-file or SAMPLE_DIR environment variable")
    if sample_id is None:
        missing_params.append("--sample-id or SAMPLE_ID environment variable")
    if batch_id is None:
        missing_params.append("--batch-id or BATCH_ID environment variable")
    if timeline_file is None:
        missing_params.append("--timeline-file or TIMELINE_FILE environment variable")
    if primer_file is None:
        missing_params.append("--primer-file or PRIMER_FILE environment variable")
    if output_fp is None:
        missing_params.append("--output-fp or RESULTS_DIR environment variable")
    if reference is None:
        missing_params.append("--reference or NEXTCLADE_REFERENCE environment variable")

    if missing_params:
        typer.echo(f"Error: Missing required parameters: {', '.join(missing_params)}")
        raise typer.Exit(1)

    typer.echo("Starting V-PIPE to SILO conversion.")

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
