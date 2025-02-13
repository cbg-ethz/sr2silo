"""Script to convert BAM to AlignedReads Objects containing all gene information.

   This script will extract nuclitide sequences and insertions from a BAM file and
    then align them to a reference genome. It will then extract amino acid sequences
    and amino acid insertions using Diamond / BlastX. Finally, it will create a JSON
    representation of the AlignedReads objects.

    That looks alike a sinle line the inputs for SILO, excluding the metadata.
"""
from __future__ import annotations

import logging
from pathlib import Path

import sr2silo.process as process
import sr2silo.process.convert as convert
import sr2silo.process.translate_align as translate_align

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def main():
    """Main function to process SAM files and generate JSON output."""
    # INPUT_NUC_ALIGMENT_FILE = Path("input/combined.bam")
    INPUT_NUC_ALIGMENT_FILE = Path("input/combined_sorted.bam")
    # INPUT_NUC_ALIGMENT_FILE = Path("input/REF_aln_trim.bam")
    FASTQ_NUC_ALIGMENT_FILE = Path("output_with_indels.fastq")
    FASTA_NUC_INSERTIONS_FILE = Path("output_ins.fasta")
    AA_ALIGNMENT_FILE = Path("diamond_blastx.sam")
    AA_REFERENCE_FILE = Path("../resources/sars-cov-2/aa_reference_genomes.fasta")
    NUC_REFERENCE_FILE = Path("../resources/sars-cov-2/nuc_reference_genomes.fasta")

    if not all(
        f.exists()
        for f in [NUC_REFERENCE_FILE, AA_REFERENCE_FILE, INPUT_NUC_ALIGMENT_FILE]
    ):
        raise FileNotFoundError("One or more input files are missing")

    # sort and index the input BAM file
    INPUT_NUC_ALIGMENT_FILE_sorted_indexed = Path("input/combined_sorted.bam")
    convert.sort_and_index_bam(
        INPUT_NUC_ALIGMENT_FILE, INPUT_NUC_ALIGMENT_FILE_sorted_indexed
    )

    logging.info("Parsing Nucliotide: BAM FASTQ conversion (with INDELS)")
    convert.bam_to_fastq_handle_indels(
        INPUT_NUC_ALIGMENT_FILE_sorted_indexed,
        FASTQ_NUC_ALIGMENT_FILE,
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

    gene_set = process.get_gene_set_from_ref(AA_REFERENCE_FILE)
    logging.info(f"Loaded gene reference with genes: {gene_set}")

    logging.info("Processing nucleotide alignments")
    aligned_reads = translate_align.read_in_AligendReads_nuc_seq(
        FASTQ_NUC_ALIGMENT_FILE, nuc_reference_length, gene_set
    )

    logging.info("Adding nucleotide insertions to reads")
    aligned_reads = translate_align.read_in_AligendReads_nuc_ins(
        aligned_reads, FASTA_NUC_INSERTIONS_FILE
    )

    # Process AA alignment file and update corresponding reads
    logging.info("Processing AA alignments")
    aligned_reads = translate_align.read_in_AlignedReads_aa_seq_and_ins(
        aligned_reads, AA_ALIGNMENT_FILE, gene_set
    )

    for read_id, read in list(aligned_reads.items()):
        print(read)


if __name__ == "__main__":
    """Run the main function."""
    main()
