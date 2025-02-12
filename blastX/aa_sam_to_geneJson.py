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
import json


from sr2silo.process.interface import AAInsertion, AAInsertionSet, AlignedRead, Gene, NucInsertion
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

    def insert_nuc_read(self, read: AlignedRead):
        """Insert a nucleotide read into the store."""
        ins_json = json.dumps([ins.__dict__ for ins in read.nucleotide_insertions])
        self.conn.execute(
            "INSERT OR REPLACE INTO reads (read_id, unaligned_seq, aligned_nuc_seq, "
            "nucleotide_insertions) VALUES (?,?,?,?)",
            (
                read.read_id,
                read.unaligned_nucleotide_sequences,
                read.aligned_nucleotide_sequences,
                ins_json,
            ),
        )
        self.conn.commit()

    def update_aa_alignment(self, read_id: str, aligned_aa_seq: str, aa_insertions: dict):
        """Update the AA alignment for a read in the store."""
        aa_ins_json = json.dumps(
            {k: [str(i) for i in v] for k, v in aa_insertions.items()}
        )
        self.conn.execute(
            "UPDATE reads SET aligned_aa_seq=?, aa_insertions=? WHERE read_id=?",
            (aligned_aa_seq, aa_ins_json, read_id),
        )
        self.conn.commit()

    def update_nuc_ins(self, record: tuple[str, NucInsertion]):
        """Update the nucleotide insertions for a read in the store."""
        read_id, nuc_ins = record
        cursor = self.conn.execute("SELECT nucleotide_insertions FROM reads WHERE read_id=?", (read_id,))
        row = cursor.fetchone()
        if row:
            nuc_ins_list = json.loads(row[0])
            nuc_ins_list.append(nuc_ins.__dict__)
            nuc_ins_json = json.dumps(nuc_ins_list)
            self.conn.execute(
                "UPDATE reads SET nucleotide_insertions=? WHERE read_id=?",
                (nuc_ins_json, read_id),
            )
            self.conn.commit()

    def bulk_insert_nuc_reads(self, records: list[AlignedRead]):
        """Bulk insert multiple nucleotide reads."""
        query = (
            "INSERT OR REPLACE INTO reads "
            "(read_id, unaligned_seq, aligned_nuc_seq, nucleotide_insertions) "
            "VALUES (?,?,?,?)"
        )
        data = [
            (
                read.read_id,
                read.unaligned_nucleotide_sequences,
                read.aligned_nucleotide_sequences,
                json.dumps([ins.__dict__ for ins in read.nucleotide_insertions]),
            )
            for read in records
        ]
        self.conn.executemany(query, data)
        self.conn.commit()

    def bulk_update_aa_alignments(self, records: list[tuple[str, AAInsertionSet, dict]]):
        """Bulk update multiple AA alignments."""
        query = "UPDATE reads SET aligned_aa_seq=?, aa_insertions=? WHERE read_id=?"
        data = [
            (
                aligned_aa_seq,
                aa_insertions.to_dict(),
                read_id,
            )
            for aligned_aa_seq, aa_insertions, read_id in records
        ]
        self.conn.executemany(query, data)
        self.conn.commit()

    # TODO: this cannot work as there are way to many reads to load into memory
    def get_all_reads(self) -> list[AlignedRead]:
        """Retrieve all reads from the store."""
        cursor = self.conn.execute("SELECT * FROM reads")
        all_reads = []
        for row in cursor:
            read = AlignedRead(
                read_id=row[0],
                unaligned_nucleotide_sequences=row[1],
                aligned_nucleotide_sequences=row[2],
                nucleotide_insertions=json.loads(row[3]) if row[3] else [],
                aligned_amino_acid_sequences=row[4],
                amino_acid_insertions=json.loads(row[5]) if row[5] else {},
            )
            all_reads.append(read)
        return all_reads

    def get_read(self, read_id: str) -> AlignedRead:
        """Retrieve a read by its ID from the store."""
        cursor = self.conn.execute("SELECT * FROM reads WHERE read_id=?", (read_id,))
        row = cursor.fetchone()
        if row:
            aa_insertions = json.loads(row[5]) if row[5] else {}
            print(aa_insertions)
            nuc_insertions = json.loads(row[3]) if row[3] else []
            print(nuc_insertions)

            return AlignedRead(
                read_id=row[0],
                unaligned_nucleotide_sequences=row[1],
                aligned_nucleotide_sequences=row[2],
                nucleotide_insertions=aa_insertions,
                aligned_amino_acid_sequences=row[4],
                amino_acid_insertions=nuc_insertions,
            )
        else:
            return None

    def dump_all_json(self) -> str:
        """Dump all reads in the store to a JSON string."""
        all_reads = self.get_all_reads()
        return json.dumps([read.__dict__ for read in all_reads], indent=4)

def nuc_alignment_to_AlignedRead(nuc_alignment_fp: Path, aa_alignment_fp: Path):

    return NotImplementedError

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
    logging.info(f"Loaded gene reference with genes: {gene_dict.keys()}")


    with tempfile.NamedTemporaryFile(delete=False) as temp_db:
        #read_store = ReadStore(db_path=temp_db.name)
        read_store = ReadStore(db_path="read_store.db")
        logging.info(f"Using temporary database at {temp_db.name}")

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
                        amino_acid_insertions = AAInsertionSet(gene_dict.items()),
                        aligned_amino_acid_sequences= None,
                    )

                    batch_nuc_records.append(read)

                    if len(batch_nuc_records) >= BATCH_SIZE:
                        read_store.bulk_insert_nuc_reads(batch_nuc_records)
                        batch_nuc_records = []
                    pbar.update(1)
                # Insert any remaining records
                if batch_nuc_records:
                    read_store.bulk_insert_nuc_reads(batch_nuc_records)

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
                print(nuc_ins_record)
                read_store.update_nuc_ins(nuc_ins_record)



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

                    batch_aa_records.append(
                        (padded_aa_alignment, aa_insertions, read_id)
                    )

                    if len(batch_aa_records) >= BATCH_SIZE:
                        read_store.bulk_update_aa_alignments(batch_aa_records)
                        batch_aa_records = []
                    pbar.update(1)
                # Update any remaining records
                if batch_aa_records:
                    read_store.bulk_update_aa_alignments(batch_aa_records)

        # Dump combined results to NDJSON incrementally to handle large data
        logging.info("Dumping final NDJSON")
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


        print("try get read")
        read =  read_store.get_read("AV233803:AV044:2411515907:1:11305:3061:3014")
        print(read.to_json())

        # read in the last JSON as AlignedRead object
        with open(final_json_fp, "r") as f:
            #for line in f:
            #    read = AlignedRead(**json.loads(line))
            #    print(read.to_json())

            # read the last line of the file
            last_line = f.readlines()[-1]
            read = AlignedRead(**json.loads(last_line))
            print(read.to_json())

        # os.remove(temp_db.name)


if __name__ == "__main__":
    """Run the main function."""
    main()
