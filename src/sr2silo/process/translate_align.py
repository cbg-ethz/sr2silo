"""Implements the translation of nucleotides alignments to amino acid alignments."""

from __future__ import annotations

import logging
import subprocess
import os
import tempfile
from pathlib import Path
from typing import List

from sr2silo.process.convert import bam_to_fasta


# TODO: to remove // perhaps to test diamond against nextclade ?
def translate_nextclade(
    input_files: List[Path], result_dir: Path, nextclade_reference: str
) -> None:
    """Translate consensus nucleotides to amino acid sequences.

    Args:
        input_file (str): The path to the input file.
                          the nucleotide sequences in fasta format.
        result_dir (str): The path to the directory to save the results.
        nextclade_reference (str): The path to the nextclade reference.
                                    e.g. nextstrain/sars-cov-2/XBB
                                    see `nextclade dataset list`
    """

    with tempfile.TemporaryDirectory() as temp_dir:
        logging.debug(f"temp_dir: {temp_dir}")
        # first get the test dataset from the gff3 file
        command = [
            "nextclade",
            "dataset",
            "get",
            "--name",
            f"{nextclade_reference}",
            "--output-dir",
            temp_dir,
        ]
        logging.debug(f"Running command: {command}")
        subprocess.run(command, check=True)

        for input_file in input_files:
            logging.info(f"Translating {input_file}")

            # then replace the sequences.fasta in the temp_dir
            #  with the sequences.fasta from the input file
            command = ["cp", input_file, f"{temp_dir}/sequences.fasta"]
            logging.debug(f"Running command: {command}")
            subprocess.run(command, check=True)

            # then run the nextclade run command
            command = [
                "nextclade",
                "run",
                "--input-dataset",
                temp_dir,
                f"--output-all={result_dir}/",
                f"{temp_dir}/sequences.fasta",
            ]

            logging.debug(f"Running nextclade: {command}")

            try:
                result = subprocess.run(
                    command, check=True, capture_output=True, text=True
                )
                logging.debug(result.stdout)
                logging.debug(result.stderr)
            except subprocess.CalledProcessError as e:
                logging.error(f"nextclade failed with exit code {e.returncode}")
                logging.error(e.stderr)
                raise

            # move the results to the result_dir
            result_path = result_dir / input_file.stem
            command = ["mv", f"{temp_dir}/results", str(result_path)]



def nuc_to_aa_alignment(
    in_nuc_alignment_fp: Path,
    in_aa_reference_fp: Path,
    out_aa_alignment_fp: Path,
) -> None:
    """
    Function to convert files and translate and align with Diamond / blastx.

    Args:
        in_nuc_alignment_fp (Path): Path to the input nucleotide alignment file.
        in_aa_reference_fp (Path): Path to the input amino acid reference file.
        out_aa_alignment_fp (Path): Path to the output amino acid alignment file.

    Returns:
        None

    Description:
        Uses Diamond with the settings:
        --evalue 1
        --gapopen 6
        --gapextend 2
        --outfmt 101
        --matrix BLOSUM62
        --unal 0
        --max-hsps 1
        --block-size 0.5
    """

    # temporary fasta file for AA alignment
    fasta_nuc_for_aa_alignment = out_aa_alignment_fp.with_suffix(".tmp.fasta")

    logging.info("Converting BAM to FASTQ for AA alignment")
    logging.info("FASTA conversion for AA alignment")
    bam_to_fasta(in_nuc_alignment_fp, fasta_nuc_for_aa_alignment)

    try:
        db_ref_fp = Path(in_aa_reference_fp.stem + ".temp.db")
        # ==== Make Sequence DB ====
        logging.info("Diamond makedb")
        print("== Making Sequence DB ==")
        result = os.system(
            f"diamond makedb --in {in_aa_reference_fp} -d {db_ref_fp}"
        )
        if result != 0:
            raise RuntimeError(
                "Error occurred while making sequence DB with diamond makedb"
            )
    except Exception as e:
        print(f"An error occurred while making sequence DB: {e}")
        raise

    try:
        # ==== Alignment ====
        logging.info("Diamond blastx alignment")
        result = os.system(
            f"diamond blastx -d {db_ref_fp} -q {fasta_nuc_for_aa_alignment} "
            f"-o {out_aa_alignment_fp} "
            f"--evalue 1 --gapopen 6 --gapextend 2 --outfmt 101 --matrix BLOSUM62 "
            f"--unal 0 --max-hsps 1 --block-size 0.5"
        )
        if result != 0:
            raise RuntimeError(
                "Error occurred while aligning to AA with diamond blastx"
            )
    except Exception as e:
        print(f"An error occurred while aligning to AA: {e}")
        raise
    finally:
        # Ensure the temporary fasta file is deleted
        if fasta_nuc_for_aa_alignment.exists():
            fasta_nuc_for_aa_alignment.unlink()

    return None