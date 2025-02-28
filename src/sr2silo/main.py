"""Entry point for the sr2silo CLI."""

from __future__ import annotations

import typer

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
def import_to_loculus():
    """
    V-PIPE to SILO conversion with amio acids, special metadata,
    Upload to S3 and submission to Loculus.
    """
    typer.echo("Starting V-PIPE to SILO conversion.")
