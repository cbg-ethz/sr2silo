"""Entry point for the sr2silo CLI."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Annotated

import typer

from sr2silo.config import is_ci_environment
from sr2silo.run import process_file

app = typer.Typer(
    name="sr2silo",
    help=(
        "Convert Short-Read nulclitide .bam alignments to cleartext alignments, "
        "with amino acids and insertions, in JSON format."
    ),
)


@app.command()
def run():
    """
    Wrangel short-reads into cleartext alignments,
    optionally translation and align in amino acids.
    """
    typer.echo("Not yet implemented.")


@app.command()
def import_to_loculus(
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
    timeline_file: Annotated[
        Path,
        typer.Option(
            "--timeline-file",
            "-t",
            help="Path to the timeline file.",
        ),
    ],
    primer_file: Annotated[
        Path,
        typer.Option(
            "--primer-file",
            "-p",
            help="Path to the primers file.",
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
    reference: Annotated[
        str,
        typer.Option(
            "--reference",
            "-r",
            help="See folder names in resources/",
        ),
    ] = "sars-cov-2",
    upload: Annotated[
        bool,
        typer.Option(
            "--upload/--no-upload",
            help="Upload and submit to SILO.",
        ),
    ] = False,
) -> None:
    """
    V-PIPE to SILO conversion with amio acids, special metadata,
    Upload to S3 and submission to Loculus.
    """
    typer.echo("Starting V-PIPE to SILO conversion.")

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
        input_file=input_file,
        sample_id=sample_id,
        batch_id=batch_id,
        timeline_file=timeline_file,
        primers_file=primer_file,
        output_fp=output_fp,
        reference=reference,
        upload=upload,
    )
