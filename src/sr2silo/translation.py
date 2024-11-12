"""This module contains the logic to translate consensus nuclotides to
    amino acid sequences."""

from __future__ import annotations

import logging
import subprocess


def translate(input_file: str, result_dir: str, nextclade_reference: str) -> None:
    """Translate consensus nucleotides to amino acid sequences.

    Args:
        input_file (str): The path to the input file.
                          the nucleotide sequences in fasta format.
        result_dir (str): The path to the directory to save the results.
        nextclade_reference (str): The path to the nextclade reference.
                                    e.g. nextstrain/sars-cov-2/XBB
                                    see `nextclade dataset list`
    """

    # first get the test dataset from the gff3 file
    command = [
        "nextclade",
        "dataset",
        "get",
        "--name",
        f"{nextclade_reference}",
        "--output-dir",
        "data/sars-cov-2",
    ]
    logging.debug(f"Running command: {command}")
    subprocess.run(command, check=True)

    # then replace the sequences.fasta in the
    # data/sars-cov-2 with the sequences.fasta from the input file
    command = ["cp", input_file, "data/sars-cov-2/sequences.fasta"]
    logging.debug(f"Running command: {command}")
    subprocess.run(command, check=True)

    # then run the nextclade run command
    command = [
        "nextclade",
        "run",
        "--input-dataset",
        "data/sars-cov-2",
        f"--output-all={result_dir}/",
        "data/sars-cov-2/sequences.fasta",
    ]

    logging.debug(f"Running nextclade: {command}")

    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        logging.debug(result.stdout)
        logging.debug(result.stderr)
    except subprocess.CalledProcessError as e:
        logging.error(f"nextclade failed with exit code {e.returncode}")
        logging.error(e.stderr)
        raise
