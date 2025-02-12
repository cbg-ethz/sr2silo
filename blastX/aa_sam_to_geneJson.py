"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations

import json
import logging
import os
import tempfile
import time
from pathlib import Path

import psutil
import pysam
from tqdm import tqdm
import json


from sr2silo.process.interface import (
    AAInsertion,
    AAInsertionSet,
    AlignedRead,
    Gene,
    NucInsertion,
    AASequenceSet,
    GeneName,
)
import sr2silo.process.convert as convert
from sr2silo.process import pad_alignment
import sr2silo.process as process

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)



def main():
    """Main function to process SAM files and generate JSON output."""
    #INPUT_NUC_ALIGMENT_FILE = "input/combined.bam"
    INPUT_NUC_ALIGMENT_FILE = "input/combined_sorted.bam"
    # INPUT_NUC_ALIGMENT_FILE = "input/REF_aln_trim.bam"
    FASTQ_NUC_ALIGMENT_FILE_WITH_INDELS = "output_with_indels.fastq"
    FASTA_NUC_INSERTIONS_FILE = "output_ins.fasta"
    AA_ALIGNMENT_FILE = "diamond_blastx.sam"
    AA_REFERENCE_FILE = "../resources/sars-cov-2/aa_reference_genomes.fasta"
    NUC_REFERENCE_FILE = "../resources/sars-cov-2/nuc_reference_genomes.fasta"
    BATCH_SIZE = 10000  # adjust as needed

    # Validate that all Resources are there before starting
    (
        NUC_REFERENCE_FILE,
        AA_REFERENCE_FILE,
        INPUT_NUC_ALIGMENT_FILE,
        AA_ALIGNMENT_FILE,
    ) = [
        Path(f)
        for f in [
            NUC_REFERENCE_FILE,
            AA_REFERENCE_FILE,
            INPUT_NUC_ALIGMENT_FILE,
            AA_ALIGNMENT_FILE
        ]
    ]
    if not all(
        f.exists()
        for f in [NUC_REFERENCE_FILE, AA_REFERENCE_FILE, INPUT_NUC_ALIGMENT_FILE]
    ):
        raise FileNotFoundError("One or more input files are missing")

    # sort and index the input BAM file
    INPUT_NUC_ALIGMENT_FILE_sorted_indexed = Path("input/combined_sorted.bam")
    convert.sort_and_index_bam(INPUT_NUC_ALIGMENT_FILE, INPUT_NUC_ALIGMENT_FILE_sorted_indexed)

    logging.info("Parsing Nucliotide: BAM FASTQ conversion (with INDELS)")
    convert.bam_to_fastq_handle_indels(
        INPUT_NUC_ALIGMENT_FILE_sorted_indexed,
        FASTQ_NUC_ALIGMENT_FILE_WITH_INDELS,
        FASTA_NUC_INSERTIONS_FILE,
    )

    # Call translation and alignment to prepare the files for downstream processing.
    process.nuc_to_aa_alignment(
        in_nuc_alignment_fp=INPUT_NUC_ALIGMENT_FILE_sorted_indexed,
        in_aa_reference_fp=AA_REFERENCE_FILE,
        out_aa_alignment_fp=AA_ALIGNMENT_FILE,
    )

    with open(NUC_REFERENCE_FILE, "r") as f:
        nuc_reference = f.read()
    nuc_reference_length = len(nuc_reference)
    logging.info(f"Loaded nucleotide reference with length {nuc_reference_length}")

    # Load gene reference
    gene_dict = convert.get_genes_and_lengths_from_ref(AA_REFERENCE_FILE)
    gene_names = [GeneName(k) for k in gene_dict.keys()]
    logging.info(f"Loaded gene reference with genes: {gene_dict.keys()}")


    # make list of AlignedRead objects
    aligned_reads = dict()

    ## Process nucleotide alignment reads incrementally
    logging.info("Processing nucleotide alignments")
    batch_nuc_records = []
    with open(FASTQ_NUC_ALIGMENT_FILE_WITH_INDELS, "r") as f:
        total_lines = sum(1 for _ in f) // 5  # Each entry consists of 5 lines
        f.seek(0)  # Reset file pointer to the beginning

        with tqdm(
            total=total_lines, desc="Processing nucleotide alignments"
        ) as pbar:
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
                    logging.error(
                        "Malformed FASTQ record encountered, skipping..."
                    )
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

                aligned_nuc_seq = pad_alignment(seq, pos, nuc_reference_length)

                read = AlignedRead(
                    read_id=read_id,
                    unaligned_nucleotide_sequences=seq,
                    aligned_nucleotide_sequences=aligned_nuc_seq,
                    nucleotide_insertions= list(),
                    amino_acid_insertions = AAInsertionSet(gene_names),
                    aligned_amino_acid_sequences= AASequenceSet(gene_names),
                )

                aligned_reads.update({read_id: read})


    ## Add Nuc insertions to the Aligned Reads incrementally
    logging.info("Adding nucleotide insertions to reads")

    with open(FASTA_NUC_INSERTIONS_FILE, "r") as f:
        # read each line seperated by tabs, read_id, position, sequence, quality
        for line in f:
            fields = line.strip().split("\t")
            read_id = fields[0]
            pos = int(fields[1])
            seq = fields[2]
            quality = fields[3] # not used

            nuc_ins = NucInsertion(position=pos, sequence=seq)
            nuc_ins_record = (read_id, nuc_ins)

            aligned_reads[read_id].nucleotide_insertions.append(nuc_ins)


    # Process AA alignment file and update corresponding reads
    logging.info("Processing AA alignments")
    batch_aa_records = []
    with open(AA_ALIGNMENT_FILE, "r") as f:
        total_lines = sum(1 for _ in f)
        f.seek(0)  # Reset file pointer to the beginning

        with tqdm(total=total_lines, desc="Processing AA alignments") as pbar:
            for line in f:
                if line.startswith("@"):  # skip header of .sam file
                    pbar.update(1)
                    continue
                fields = line.strip().split("\t")
                read_id = fields[0]
                gene_name = fields[2]
                pos = int(fields[3])
                cigar = fields[5]
                seq = fields[9]
                (
                    aa_aligned,
                    aa_insertions,
                    aa_deletions,
                ) = convert.sam_to_seq_and_indels(seq, cigar)

                padded_aa_alignment = pad_alignment(
                    aa_aligned, pos, gene_dict[gene_name].gene_length
                )

                ## update the insertions set with the new insertions

                aa_ins = AAInsertion(position=pos, sequence=aa_insertions)

                aligned_reads[read_id].amino_acid_insertions.set_insertions_for_gene(gene_name, aa_ins)
                aligned_reads[read_id].aligned_amino_acid_sequences.set_sequence(gene_name, padded_aa_alignment)

    for read_id, read in list(aligned_reads.items())[-3:]:
        print(read.to_json())


if __name__ == "__main__":
    """Run the main function."""
    main()
