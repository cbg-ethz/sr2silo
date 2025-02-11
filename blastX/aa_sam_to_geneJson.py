"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations

import json
import logging
import os
import sqlite3
import tempfile
import time
from pathlib import Path

import psutil
import pysam
from tqdm import tqdm

import sr2silo.process.convert as convert
from sr2silo.process import pad_alignment
import sr2silo.process as process

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


class ReadStore:
    """Class to store and manage reads sequences and indels in a SQLite database."""

    def __init__(self, db_path="read_store.db"):
        """Initialize the ReadStore with a SQLite database."""
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
            "INSERT OR REPLACE INTO reads (read_id, unaligned_seq, aligned_nuc_seq, "
            "nucleotide_insertions) VALUES (?,?,?,?)",
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

    def bulk_insert_nuc_reads(self, records):
        """Bulk insert multiple nucleotide reads."""
        query = (
            "INSERT OR REPLACE INTO reads "
            "(read_id, unaligned_seq, aligned_nuc_seq, nucleotide_insertions) "
            "VALUES (?,?,?,?)"
        )
        self.conn.executemany(query, records)
        self.conn.commit()

    def bulk_update_aa_alignments(self, records):
        """Bulk update multiple AA alignments."""
        query = "UPDATE reads SET aligned_aa_seq=?, aa_insertions=? WHERE read_id=?"
        self.conn.executemany(query, records)
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
        """Initialize the PerfMonitor with a label."""
        self.label = label

    def __enter__(self):
        logging.info(f"Starting {self.label}")
        self.start_time = time.time()
        self.start_mem = psutil.Process(os.getpid()).memory_info().rss
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        end_time = time.time()
        end_mem = psutil.Process(os.getpid()).memory_info().rss
        mem_diff_mb = (end_mem - self.start_mem) / (1024 * 1024)
        elapsed = end_time - self.start_time

        logging.info(
            f"{self.label} completed in {elapsed:.2f}s "
            f"with memory change {mem_diff_mb:.2f} MB"
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

    with PerfMonitor("Parsing Nucliotide: BAM FASTQ conversion (with INDELS)"):
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

    gene_dict = convert.get_genes_and_lengths_from_ref(AA_REFERENCE_FILE)
    logging.info(f"Loaded gene reference with genes: {gene_dict.keys()}")

    with tempfile.NamedTemporaryFile(delete=False) as temp_db:
        read_store = ReadStore(db_path=temp_db.name)

        ## Process nucleotide alignment reads incrementally
        with PerfMonitor("Processing nucleotide alignments"):
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
                        ins_json = json.dumps([])  # empty nucleotide insertions list
                        batch_nuc_records.append(
                            (read_id, seq, aligned_nuc_seq, ins_json)
                        )

                        if len(batch_nuc_records) >= BATCH_SIZE:
                            read_store.bulk_insert_nuc_reads(batch_nuc_records)
                            batch_nuc_records = []
                        pbar.update(1)
                    # Insert any remaining records
                    if batch_nuc_records:
                        read_store.bulk_insert_nuc_reads(batch_nuc_records)

        # Process AA alignment file and update corresponding reads
        with PerfMonitor("Processing AA alignments"):
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
                        # Build a dict for AA insertions with all genes as keys
                        aa_ins_dict = {gene: [] for gene in gene_dict.keys()}
                        aa_ins_dict[gene_name] = aa_insertions
                        padded_aa_alignment = pad_alignment(
                            aa_aligned, pos, gene_dict[gene_name].gene_length
                        )
                        aa_ins_json = json.dumps(
                            {k: [str(i) for i in v] for k, v in aa_ins_dict.items()}
                        )
                        batch_aa_records.append(
                            (padded_aa_alignment, aa_ins_json, read_id)
                        )

                        if len(batch_aa_records) >= BATCH_SIZE:
                            read_store.bulk_update_aa_alignments(batch_aa_records)
                            batch_aa_records = []
                        pbar.update(1)
                    # Update any remaining records
                    if batch_aa_records:
                        read_store.bulk_update_aa_alignments(batch_aa_records)

        # Dump combined results to NDJSON incrementally to handle large data
        with PerfMonitor("Dumping final NDJSON"):
            final_json_fp = "reads.ndjson"

            with open(final_json_fp, "w") as f:
                cursor = read_store.conn.execute("SELECT * FROM reads")
                for row in cursor:
                    read_obj = {
                    "read_id": row[0],
                    "unaligned_nucleotide_sequences": row[1],
                    "aligned_nucleotide_sequences": row[2],
                    "nucleotide_insertions": json.loads(row[3]) if row[3] else [],
                    "aligned_amino_acid_sequences": row[4],
                    "amino_acid_insertions": json.loads(row[5]) if row[5] else {},
                    }
                    f.write(json.dumps(read_obj) + "\n")

        # print the final json last element
        with open(final_json_fp, "r") as f:
            lines = f.readlines()
            if lines:
                print(lines[-1])

        print("Done!")

        os.remove(temp_db.name)


if __name__ == "__main__":
    """Run the main function."""
    main()
