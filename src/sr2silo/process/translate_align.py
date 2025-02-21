"""Implements the translation of nucleotides alignments to amino acid alignments."""

from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List

from tqdm import tqdm

import sr2silo.process.convert as convert
from sr2silo.process.interface import (
    AAInsertionSet,
    AASequenceSet,
    AlignedRead,
    GeneName,
    GeneSet,
    NucInsertion,
)


# TODO: to use as orthogonal test against blastX
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
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)

        # temporary fasta file for AA alignment
        fasta_nuc_for_aa_alignment = temp_dir_path / out_aa_alignment_fp.with_suffix(
            ".fasta"
        )

        logging.info("Converting BAM to FASTQ for AA alignment")
        logging.info("FASTA conversion for AA alignment")
        convert.bam_to_fasta(in_nuc_alignment_fp, fasta_nuc_for_aa_alignment)

        # temporary file file for amino acid reference DB
        db_ref_fp = temp_dir_path / Path(in_aa_reference_fp.stem + ".temp.db")
        try:
            # ==== Make Sequence DB ====
            logging.info("Diamond makedb")
            logging.info("== Making Sequence DB ==")
            result = os.system(
                f"diamond makedb --in {in_aa_reference_fp} -d {db_ref_fp}"
            )
            if result != 0:
                raise RuntimeError(
                    f"Error occurred while making sequence DB with diamond makedb "
                    f"- Error Code: {result}"
                )
        except Exception as e:
            logging.error(
                f"An error occurred while making sequence DB - Error Code: {e}"
            )
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
            logging.error(f"An error occurred while aligning to AA: {e}")
            raise
        finally:
            # Ensure the temporary fasta file is deleted
            if fasta_nuc_for_aa_alignment.exists():
                fasta_nuc_for_aa_alignment.unlink()

    return None


def read_in_AlignedReads_nuc_seq(
    fastq_nuc_alignment_file: Path, nuc_reference_length: int, gene_set: GeneSet
) -> Dict[str, AlignedRead]:
    """Read aligned reads from a FASTQ file with indels.

    Args:
        fastq_nuc_alignment_file (Path): Path to the FASTQ file with alignment
                                        positions, produced by
                                        `bam_to_fastq_handle_indels`.
        nuc_reference_length (int): Length of the nucleotide reference genome.
        gene_set (GeneSet): Set of genes for amino acid sequence alignment.

    Returns:
        dict[str, AlignedRead]: Dictionary of read IDs to AlignedRead objects.
    """
    aligned_reads = dict()
    with open(fastq_nuc_alignment_file, "r") as f:
        total_lines = sum(1 for _ in f) // 5  # Each entry consists of 5 lines
        f.seek(0)  # Reset file pointer to the beginning
        with tqdm(total=total_lines, desc="Processing nucleotide alignments"):
            while True:
                lines = [f.readline().strip() for _ in range(5)]
                if not lines[0]:
                    break  # End of file
                if not lines[0].startswith("@"):
                    continue
                if not (
                    lines[0].startswith("@")
                    and lines[2].startswith("+")
                    and lines[4].startswith("alignment_position:")
                ):
                    logging.error("Malformed FASTQ record encountered, skipping...")
                    continue
                read_id = lines[0][1:]
                seq = lines[1]
                try:
                    pos = int(lines[4].split(":", 1)[1])
                except Exception as e:
                    logging.error(
                        f"Error parsing alignment position for {read_id}: {e}"
                    )
                    continue
                aligned_nuc_seq = convert.pad_alignment(seq, pos, nuc_reference_length)
                read = AlignedRead(
                    read_id=read_id,
                    unaligned_nucleotide_sequences=seq,
                    aligned_nucleotide_sequences=aligned_nuc_seq,
                    nucleotide_insertions=list(),
                    amino_acid_insertions=AAInsertionSet(gene_set.get_gene_name_list()),
                    aligned_amino_acid_sequences=AASequenceSet(
                        gene_set.get_gene_name_list()
                    ),
                )
                aligned_reads.update({read_id: read})
    return aligned_reads


def read_in_AlignedReads_nuc_ins(
    aligned_reads: dict[str, AlignedRead], fasta_nuc_insertions_file: Path
) -> Dict[str, AlignedRead]:
    """Read in nucleotide insertions from a FASTA file and update the reads."""

    with open(fasta_nuc_insertions_file, "r") as f:
        # read each line separated by tabs, read_id, position, sequence, quality
        for line in f:
            fields = line.strip().split("\t")
            read_id = fields[0]
            pos = int(fields[1])
            seq = fields[2]
            # quality = fields[3]
            nuc_ins = NucInsertion(position=pos, sequence=seq)

            aligned_reads[read_id].set_nuc_insertion(nuc_ins)

    return aligned_reads


def read_in_AlignedReads_aa_seq_and_ins(
    aligned_reads: Dict[str, AlignedRead],
    fasta_aa_alignment_file: Path,
    gene_set: GeneSet,
) -> Dict[str, AlignedRead]:
    """Read in amino acid sequences and insertions from a FASTA file"""
    with open(fasta_aa_alignment_file, "r") as f:
        total_lines = sum(1 for _ in f)
        f.seek(0)  # Reset file pointer to the beginning
        with tqdm(total=total_lines, desc="Processing AA alignments") as pbar:
            for line in f:
                if line.startswith("@"):  # skip header of .sam file
                    pbar.update(1)
                    continue
                fields = line.strip().split("\t")
                read_id = fields[0]
                gene_name = GeneName(fields[2])
                pos = int(fields[3])
                cigar = fields[5]
                seq = fields[9]
                (
                    aa_aligned,
                    aa_insertions,
                    aa_deletions,
                ) = convert.sam_to_seq_and_indels(seq, cigar)
                padded_aa_alignment = convert.pad_alignment(
                    aa_aligned, pos, gene_set.get_gene_length(gene_name)
                )

                aligned_reads[read_id].amino_acid_insertions.set_insertions_for_gene(
                    gene_name, aa_insertions
                )
                aligned_reads[read_id].aligned_amino_acid_sequences.set_sequence(
                    gene_name, padded_aa_alignment
                )
    return aligned_reads


def parse_translate_align(
    nuc_reference_fp: Path, aa_reference_fp: Path, nuc_alignment_fp: Path
) -> Dict[str, AlignedRead]:
    """Parse nucleotides, translate and align amino acids the input files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)
        BAM_NUC_ALIGNMENT_FILE = temp_dir_path / "combined_sorted.bam"
        FASTQ_NUC_ALIGNMENT_FILE = temp_dir_path / "output_with_indels.fastq"
        FASTA_NUC_INSERTIONS_FILE = temp_dir_path / "output_ins.fasta"
        AA_ALIGNMENT_FILE = temp_dir_path / "diamond_blastx.sam"

        missing_files = [
            str(f)
            for f in [nuc_reference_fp, aa_reference_fp, nuc_alignment_fp]
            if not f.exists()
        ]
        if missing_files:
            raise FileNotFoundError(f"Missing input files: {', '.join(missing_files)}")

        convert.sort_and_index_bam(nuc_alignment_fp, BAM_NUC_ALIGNMENT_FILE)

        logging.info("Parsing Nucleotides: BAM FASTQ conversion (with INDELS)")
        convert.bam_to_fastq_handle_indels(
            bam_file=BAM_NUC_ALIGNMENT_FILE,
            out_fastq_fp=FASTQ_NUC_ALIGNMENT_FILE,
            out_insertions_fp=FASTA_NUC_INSERTIONS_FILE,
        )

        nuc_to_aa_alignment(
            in_nuc_alignment_fp=BAM_NUC_ALIGNMENT_FILE,
            in_aa_reference_fp=aa_reference_fp,
            out_aa_alignment_fp=AA_ALIGNMENT_FILE,
        )

        with open(nuc_reference_fp, "r") as f:
            nuc_reference = f.read()
        nuc_reference_length = len(nuc_reference)
        logging.info(f"Loaded nucleotide reference with length {nuc_reference_length}")

        gene_set = convert.get_gene_set_from_ref(aa_reference_fp)
        logging.info(f"Loaded gene reference with genes: {gene_set}")

        logging.info("Processing nucleotide alignments")
        aligned_reads = read_in_AlignedReads_nuc_seq(
            FASTQ_NUC_ALIGNMENT_FILE, nuc_reference_length, gene_set
        )

        logging.info("Adding nucleotide insertions to reads")
        aligned_reads = read_in_AlignedReads_nuc_ins(
            aligned_reads, FASTA_NUC_INSERTIONS_FILE
        )

        # Process AA alignment file and update corresponding reads
        logging.info("Processing AA alignments")
        aligned_reads = read_in_AlignedReads_aa_seq_and_ins(
            aligned_reads, AA_ALIGNMENT_FILE, gene_set
        )

    return aligned_reads


def enrich_AlignedReads_with_metadata(
    aligned_reads: Dict[str, AlignedRead],
    metadata_fp: Path,
) -> Dict[str, AlignedRead]:
    """Enrich the AlignedReads with metadata from a TSV file."""

    try:
        with open(metadata_fp, "r") as file:
            metadata = json.load(file)
    except FileNotFoundError:
        logging.error("Error: File not found")
        raise FileNotFoundError
    except json.JSONDecodeError as e:
        logging.error("Error: Invalid JSON format")
        raise json.JSONDecodeError(e.msg, e.doc, e.pos)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")
        raise e

    if metadata:
        for read_id, read in aligned_reads.items():
            read.metadata = metadata
    else:
        logging.error("No metadata found in the file")
        raise ValueError("No metadata found in the file")
    return aligned_reads


def parse_translate_align_in_batches(
    fastq_nuc_alignment_file: Path,
    nuc_reference_length: int,
    gene_set: GeneSet,
    batch_size: int = 10000,
) -> Dict[str, AlignedRead]:
    """Parse nucleotides, translate and align amino acids in batches.

    Args:
        fastq_nuc_alignment_file (Path): Path to the FASTQ file with alignment
                                        positions, produced by
                                        `bam_to_fastq_handle_indels`.
        nuc_reference_length (int): Length of the nucleotide reference genome.
        gene_set (GeneSet): Set of genes for amino acid sequence alignment.
        batch_size (int): Number of reads to process in each batch.
    """

    # split the input file into batches

    raise NotImplementedError("This function is not yet implemented")
