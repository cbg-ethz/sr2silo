"""Interactions with the Lapis API."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import requests


class LapisClient:
    """Client for interacting with the Lapis API."""

    def __init__(self, lapisUrl) -> None:
        """Initialize the Lapis client.

        Args:
            lapisUrl: Base URL for the Lapis API
        """
        self.lapisUrl = lapisUrl

    def referenceGenome(self) -> dict:
        """Fetch reference genome from the Lapis `sample/referenceGenome` endpoint.

        Returns:
            JSON response as a dictionary.

        Raises:
            Exception: If the request fails.
        """
        url = f"{self.lapisUrl}/sample/referenceGenome?downloadAsFile=false"
        headers = {"accept": "application/json"}

        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            logging.error(f"Error fetching reference genome: {e}")
            raise Exception(f"Failed to fetch reference genome: {e}")

    @staticmethod
    def referenceGenomeToFasta(
        reference_json_string: str, nucleotide_out_fp: Path, amino_acid_out_fp: Path
    ) -> None:
        """Convert a reference JSON from `sample/referenceGenome` endpoint
            to separate nucleotide and amino acid reference FASTA files.

        Args:
            reference_json_string: JSON string containing reference sequences with
                                    'nucleotideSequences' and 'genes' sections
            nucleotide_out_fp: Path to the output nucleotide FASTA file
            amino_acid_out_fp: Path to the output amino acid FASTA file

        Returns:
            None
        """
        # Parse the JSON string
        reference_data = json.loads(reference_json_string)

        # Create nucleotide FASTA file
        with nucleotide_out_fp.open("w") as nuc_file:
            for nuc_seq in reference_data.get("nucleotideSequences", []):
                seq_name = nuc_seq.get("name", "main")
                sequence = nuc_seq.get("sequence", "")
                nuc_file.write(f">{seq_name}\n{sequence}\n")

        logging.info(f"Nucleotide FASTA file created at: {nucleotide_out_fp}")

        # Create amino acid FASTA file
        with amino_acid_out_fp.open("w") as aa_file:
            for gene in reference_data.get("genes", []):
                gene_name = gene.get("name", "")
                sequence = gene.get("sequence", "")
                # Remove stop codon asterisk if present
                sequence = sequence.rstrip("*")
                aa_file.write(f">{gene_name}\n{sequence}\n")

        logging.info(f"Amino acid FASTA file created at: {amino_acid_out_fp}")
