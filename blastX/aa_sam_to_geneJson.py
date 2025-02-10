"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations

import json
import logging
import os
import sqlite3
import time
from pathlib import Path
from typing import Dict
import tempfile
import json
from tqdm import tqdm

import psutil
from sr2silo.process import pad_alignment
import sr2silo.process.convert as convert


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


class ReadStore:
    """Class to store and manage reads sequences and indels in a SQLite database."""
    def __init__(self, db_path="read_store.db"):
        self.conn = sqlite3.connect(db_path)
        self.conn.execute(
            "CREATE TABLE IF NOT EXISTS reads ("
            "read_id TEXT PRIMARY KEY, "
            "unaligned_seq TEXT, "
            "aligned_nuc_seq TEXT, "
            "nucleotide_insertions TEXT, "
            "aligned_aa_seq TEXT, "
            "aa_insertions TEXT)"
        )

    def insert_nuc_read(self, read_id, unaligned_seq, aligned_nuc_seq, nuc_ins):
        """Insert a nucleotide read into the store."""
        import json

        ins_json = json.dumps([ins.__dict__ for ins in nuc_ins])
        self.conn.execute(
            "INSERT OR REPLACE INTO reads (read_id, unaligned_seq, aligned_nuc_seq, nucleotide_insertions) VALUES (?,?,?,?)",
            (read_id, unaligned_seq, aligned_nuc_seq, ins_json),
        )
        self.conn.commit()

    def update_aa_alignment(self, read_id, aligned_aa_seq, aa_insertions):
        """Update the AA alignment for a read in the store."""
        import json

        aa_ins_json = json.dumps(
            {k: [str(i) for i in v] for k, v in aa_insertions.items()}
        )
        self.conn.execute(
            "UPDATE reads SET aligned_aa_seq=?, aa_insertions=? WHERE read_id=?",
            (aligned_aa_seq, aa_ins_json, read_id),
        )
        self.conn.commit()

    def dump_all_json(self):
        """Dump all reads in the store to a JSON string."""
        cursor = self.conn.execute("SELECT * FROM reads")
        all_reads = []
        for row in cursor:
            read_obj = {
                "read_id": row[0],
                "unaligned_nucleotide_sequences": row[1],
                "aligned_nucleotide_sequences": row[2],
                "nucleotide_insertions": json.loads(row[3]) if row[3] else [],
                "aligned_amino_acid_sequences": row[4],
                "amino_acid_insertions": json.loads(row[5]) if row[5] else {},
            }
            all_reads.append(read_obj)
        return json.dumps(all_reads, indent=4)


class PerfMonitor:
    """Context manager that logs runtime and memory usage for a code block."""

    def __init__(self, label: str):
        self.label = label

    def __enter__(self):
        logging.info(f"Starting {self.label}")
        self.start_time = time.time()
        self.start_mem = psutil.Process(os.getpid()).memory_info().rss
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.time()
        end_mem = psutil.Process(os.getpid()).memory_info().rss
        elapsed = end_time - self.start_time
        mem_diff_mb = (end_mem - self.start_mem) / (1024 * 1024)
        logging.info(
            f"{self.label} completed in {elapsed:.2f}s with memory change {mem_diff_mb:.2f} MB"
        )


####################################################################################################
# Main function


def main():
    INPUT_NUC_ALIGMENT_FILE = "input/combined.bam"
    #INPUT_NUC_ALIGMENT_FILE = "input/REF_aln.bam"
    FASTA_NUC_FOR_AA_ALINGMENT = "output.fasta"
    FASTQ_NUC_ALIGMENT_FILE_WITH_INDELS = "output_with_indels.fastq"
    FASTA_NUC_INSERTIONS_FILE = "output_ins.fasta"
    AA_ALIGNMENT_FILE = "diamond_blastx.sam"
    AA_REFERENCE_FILE = "../resources/sars-cov-2/aa_reference_genomes.fasta"
    NUC_REFERENCE_FILE = "../resources/sars-cov-2/nuc_reference_genomes.fasta"

    # Validate that all Resouces are there before starting
    NUC_REFERENCE_FILE, AA_REFERENCE_FILE, INPUT_NUC_ALIGMENT_FILE = [
        Path(f)
        for f in [NUC_REFERENCE_FILE, AA_REFERENCE_FILE, INPUT_NUC_ALIGMENT_FILE]
    ]
    if not all(
        f.exists()
        for f in [NUC_REFERENCE_FILE, AA_REFERENCE_FILE, INPUT_NUC_ALIGMENT_FILE]
    ):
        raise FileNotFoundError("One or more input files are missing")

    # Sort the BAM file
    with PerfMonitor("BAM sorting"):
        convert.sort_bam_file(INPUT_NUC_ALIGMENT_FILE, Path("input/sorted.bam"))

    # Create index for BAM file
    with PerfMonitor("BAM indexing"):
        bam_file = Path("input/sorted.bam")
        bai_file = bam_file.with_suffix(".bai")
        if not bai_file.exists() or bam_file.stat().st_mtime > bai_file.stat().st_mtime:
            print("Creating index for BAM file")
            convert.create_index(bam_file)

    print("Converting BAM to FASTQ with INDELS to show NUC alignment")

    # make a mofidifed path for indels
    with PerfMonitor("FASTQ conversion (with INDELS)"):
        convert.bam_to_fastq_handle_indels(
            "input/sorted.bam",
            FASTQ_NUC_ALIGMENT_FILE_WITH_INDELS,
            FASTA_NUC_INSERTIONS_FILE,
        )

    print("Converting BAM to FASTQ for AA alignment")
    with PerfMonitor("FASTA conversion for AA alignment"):
        convert.bam_to_fasta("input/sorted.bam", FASTA_NUC_FOR_AA_ALINGMENT)

    try:
        # translate and align to AA
        # ==== Make Sequence DB ====
        with PerfMonitor("Diamond makedb"):
            print("== Making Sequence DB ==")
            result = os.system(
                f"diamond makedb --in {AA_REFERENCE_FILE} -d ref/hxb_pol_db"
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
        with PerfMonitor("Diamond blastx alignment"):
            print("== Aligning to AA ==")
            # Replace batch splitting with a single diamond blastx call on the entire FASTA file.
            result = os.system(
                f"diamond blastx -d ref/hxb_pol_db -q {FASTA_NUC_FOR_AA_ALINGMENT} -o {AA_ALIGNMENT_FILE} "
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

    with open(NUC_REFERENCE_FILE, "r") as f:
        nuc_reference = f.read()
    nuc_reference_length = len(nuc_reference)
    logging.info(f"Loaded nucleotide reference with length {nuc_reference_length}")

    gene_dict = convert.get_genes_and_lengths_from_ref(AA_REFERENCE_FILE)
    logging.info(f"Loaded gene refernece with genes: {gene_dict.keys()}")

    with tempfile.NamedTemporaryFile(delete=False) as temp_db:
        read_store = ReadStore(db_path=temp_db.name)

        ## Process nucleotide alignment reads incrementally
        with PerfMonitor("Processing nucleotide alignments"):
            with open(FASTQ_NUC_ALIGMENT_FILE_WITH_INDELS, "r") as f:
                total_lines = sum(1 for _ in f) // 5  # Each entry consists of 5 lines
                f.seek(0)  # Reset file pointer to the beginning

                with tqdm(total=total_lines, desc="Processing nucleotide alignments") as pbar:
                    while True:
                        header = f.readline().strip()
                        if not header:
                            break  # End of file
                        # check it starts with @, if not skip
                        if not header.startswith("@"):
                            continue
                        seq = f.readline().strip()
                        plus = f.readline().strip()
                        # NB @ future: propagate the quality scores from here
                        f.readline().strip()
                        pos_line = f.readline().strip()

                        if not (
                            header.startswith("@")
                            and plus.startswith("+")
                            and pos_line.startswith("alignment_position:")
                        ):
                            logging.error("Malformed FASTQ record encountered, skipping...")
                            continue

                        read_id = header[1:]
                        try:
                            pos = int(pos_line.split(":", 1)[1])
                        except Exception as e:
                            logging.error(
                            f"Error parsing alignment position for {read_id}: {e}"
                            )
                            continue

                        aligned_nuc_seq = pad_alignment(seq, pos, nuc_reference_length)
                        nuc_ins = (
                            []
                        )  # Here, keep an empty list for nucleotide insertions if needed
                        read_store.insert_nuc_read(read_id, seq, aligned_nuc_seq, nuc_ins)
                        pbar.update(1)  # Update progress bar for each processed entry

        # Process AA alignment file and update corresponding reads
        with PerfMonitor("Processing AA alignments"):
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
                        aa_aligned, aa_insertions, aa_deletions = convert.sam_to_seq_and_indels(
                            seq, cigar
                        )
                        # Build a dict for AA insertions with empty entries for other genes
                        aa_ins_dict = {gene: [] for gene in gene_dict.keys()}
                        aa_ins_dict[gene_name] = aa_insertions
                        padded_aa_alignment = pad_alignment(
                            aa_aligned, pos, gene_dict[gene_name].gene_length
                        )
                        # Update the read record with AA alignment info.
                        read_store.update_aa_alignment(
                            read_id, padded_aa_alignment, aa_ins_dict
                        )
                        pbar.update(1)  # Update progress bar for each processed line

        # Dump combined results to JSON (or write out incrementally)
        with PerfMonitor("Dumping final JSON"):
            final_json = read_store.dump_all_json()

            final_json_fp = "reads.json"

            with open(final_json_fp, "w") as f:
                f.write(final_json)

        # print the first and last AlignedRead objects
        print("First read:")
        print(json.dumps(json.loads(final_json)[0], indent=4))
        print("Last read:")
        print(json.dumps(json.loads(final_json)[-1], indent=4))

        print("Done!")

        os.remove(temp_db.name)



if __name__ == "__main__":
    main()
