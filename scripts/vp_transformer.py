"""Converts V-PIPE's outputs of nucleotide read sequencs output to
   ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import csv
import datetime
import json
import logging
from pathlib import Path

import click
import yaml

import silo_input_transformer
from sr2silo.convert import bam_to_sam
from sr2silo.lapis import submit
from sr2silo.process import pair_normalize_reads
from sr2silo.s3 import compress_bz2, upload_file_to_s3
from sr2silo.translation import translate

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def load_config(config_file: Path) -> dict:
    """Load a JSON configuration file."""
    try:
        with config_file.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"Config file not found: {config_file}")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from config file: {config_file} - {e}")
        raise


def sample_id_decoder(sample_id: str) -> dict:
    """Decode the sample ID into individual components.

    Args:
        sample_id (str): The sample ID to decode.

    Returns:
        dict: A dictionary containing the decoded components.
              containing the following keys:
                - sequencing_well_position (str : sequencing well position)
                - location_code (int : code of the location)
                - sampling_date (str : date of the sampling)
    """
    components = sample_id.split("_")
    # Assign components to meaningful variable names
    well_position = components[0]  # A1
    location_code = components[1]  # 10
    sampling_date = f"{components[2]}-{components[3]}-{components[4]}"  # 2024-09-30
    return {
        "sequencing_well_position": well_position,
        "location_code": location_code,
        "sampling_date": sampling_date,
    }


def batch_id_decoder(batch_id: str) -> dict:
    """Decode the batch ID into individual components.

    Args:
        batch_id (str): The batch ID to decode.

    Returns:
        dict: A dictionary containing the decoded components.
              containing the following keys:
                - sequencing_date (str : date of the sequencing)
                - flow_cell_serial_number (str : serial number of the flow cell)
    """
    components = batch_id.split("_")
    # Assign components to meaningful variable names
    sequencing_date = (
        f"{components[0][:4]}-{components[0][4:6]}-{components[0][6:]}"  # 2024-10-18
    )
    flow_cell_serial_number = components[1]  # AAG55WNM5
    return {
        "sequencing_date": sequencing_date,
        "flow_cell_serial_number": flow_cell_serial_number,
    }


def convert_to_iso_date(date: str) -> str:
    """Convert a date string to ISO 8601 format (date only)."""
    # Parse the date string
    date_obj = datetime.datetime.strptime(date, "%Y-%m-%d")
    # Format the date as ISO 8601 (date only)
    return date_obj.date().isoformat()


def get_metadata(sample_id: str, batch_id: str, timeline: Path, primers: Path) -> dict:
    """
    Get metadata for a given sample and batch directory.
    Cross-references the directory with the timeline file to get the metadata.

    Args:
        sample_id (str): The sample ID to use for metadata.
        batch_id (str): The batch ID to use for metadata.
        timeline (Path): The timeline file to cross-reference the metadata.
        primers (Path): The primers file to cross-reference the metadata.

    Returns:
        dict: A dictionary containing the metadata.

    """

    metadata = {}
    metadata["sample_id"] = sample_id
    metadata["batch_id"] = batch_id

    # Decompose the ids into individual components
    logging.info(f"Decoding sample_id: {metadata['sample_id']}")
    sample_id = metadata["sample_id"]
    metadata.update(sample_id_decoder(sample_id))
    logging.info(f"Decoding batch_id: {metadata['batch_id']}")
    batch_id = metadata["batch_id"]
    metadata.update(batch_id_decoder(batch_id))

    # Read the timeline file to get additional metadata
    # find row with matching sample_id and batch_id
    # timline has headers:
    #  sample_id	batch_id	read_length	primer_protocol	location_code	sampling_date	location_name
    # get read length, primer protocol, location name
    # double check if location code and location code are the same
    if not timeline.is_file():
        logging.error(f"Timeline file not found or is not a file: {timeline}")
        raise FileNotFoundError(f"Timeline file not found or is not a file: {timeline}")
    with timeline.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0] == metadata["sample_id"] and row[1] == metadata["batch_id"]:
                logging.info(
                    f"Enriching metadata with timeline data e.g. read_length, primer_protocol, location_name"
                )
                metadata["read_length"] = row[2]
                metadata["primer_protocol"] = row[3]
                metadata["location_name"] = row[6]
                # Convert sampling_date to ISO format for comparison
                timeline_sampling_date = convert_to_iso_date(row[5])
                if int(metadata["location_code"]) != int(row[4]):
                    # output both location codes for comparison and their types for debugging
                    logging.warning(
                        f"Mismatch in location code for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
                    )
                    logging.debug(
                        f"Location code mismatch: {metadata['location_code']} (sample_id) vs {row[4]} (timeline)"
                    )
                    logging.debug(
                        f"Location code types: {type(metadata['location_code'])} (sample_id) vs {type(row[4])} (timeline)"
                    )
                if metadata["sampling_date"] != timeline_sampling_date:
                    # output both sampling dates for comparison and their types for debugging
                    logging.warning(
                        f"Mismatch in sampling date for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
                    )
                    logging.debug(
                        f"Sampling date mismatch: {metadata['sampling_date']} (sample_id) vs {timeline_sampling_date} (timeline)"
                    )
                    logging.debug(
                        f"Sampling date types: {type(metadata['sampling_date'])} (sample_id) vs {type(timeline_sampling_date)} (timeline)"
                    )
                break
        else:
            raise ValueError(
                f"No matching entry found in timeline for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
            )
    # Read the primers yaml to get additional metadata
    # find the key with matching primer_protocol and get the "name" value
    # as the canonical name of the primer protocol
    if not primers.is_file():
        logging.error(f"Primers file not found or is not a file: {primers}")
        raise FileNotFoundError(f"Primers file not found or is not a file: {primers}")
    # Load YAML file
    with open(primers, "r") as file:
        primers = yaml.safe_load(file)
    logging.debug(f"Primers: {primers}")
    logging.debug(f" Type of primers: {type(primers)}")
    for primer in primers.keys():
        if primer == metadata["primer_protocol"]:
            logging.info(
                f"Enriching metadata with primer data e.g. primer_protocol_name"
            )
            metadata["primer_protocol_name"] = primers[primer]["name"]
            break
    else:
        raise ValueError(
            f"No matching entry found in primers for primer_protocol {metadata['primer_protocol']}"
        )
    return metadata


def wrangle_for_transformer(
    input_dir: Path,
    output_dir: Path,
    fasta_file: Path,
    insertions_file: Path,
    metadata_file: Path,
    database_config: Path,
) -> dict[str, Path]:
    """Wrangle the sequences to the format required by the silo_input_transformer.

    Args:
        input_dir (Path): The directory containing the sequences to wrangle.
                                aa_insertions.tsv, gene_*.fasta, metadata.json,
                                nuc_main.fasta, nucleotide_insertions.tsv,
                                reference_genomes.json, unaligned_main.fasta
        output_dir (Path): The empty working and output directory to save
                            intermediate files and results i.e. `ndjson` file.
        fasta_file (Path): The FASTA file containing the unaligned nuclotide
                            sequences.
        insertions_file (Path): The tsv file containing the nucleotide insertions.
        metadata_file (Path): The metadata json file containing the per sequencing run metadata,
                              that is copied for each read_id in the output.
        database_config (Path): The database configuration file containing the
                                schema for the metadata, has to match the
                                metadata keys in the metadata.json file but the
                                read_id key is not in the schema.

    Returns:
        dict[str, Path]: A dictionary containing the paths to the files
                         created during the wrangling process.

                         metadata_fp: metadata.tsv
                         database_config_fp: database_config.yaml
                         reference_genomes_fp: reference_genomes.json
    """

    logging.info(f"Wrangling sequences to Nextclade format")
    # make a directory for the nextclade-like results
    output_dir.mkdir(parents=True, exist_ok=True)
    # copy over the nucleotide sequnecefile and name nuc_main.fasta
    nuc_main = output_dir / "nuc_main.fasta"
    with fasta_file.open() as f:
        nuc_main_data = f.read()
    with nuc_main.open("w") as f:
        f.write(nuc_main_data)
    # copy over the amino acids sequences and name them gene_ and keep the file ending and last part of the name
    # e.g. nextclade.cds_translation.ORF1a.fasta -> gene_ORF1a.fasta
    gene_names = []
    for file in input_dir.glob("nextclade.cds_translation.*.fasta"):
        gene_name = file.name.split(".")[2]
        gene_names.append(gene_name)
        gene_file = output_dir / f"gene_{gene_name}.fasta"
        with file.open() as f:
            gene_data = f.read()
        with gene_file.open("w") as f:
            f.write(gene_data)
    # copy over the nucoletide insertions file and name it nucolotide_insertions.tsv
    # add a header of "read_id  | main" to the file
    nuc_insertions = output_dir / "nucleotide_insertions.tsv"
    with insertions_file.open() as f:
        insertions = f.read()
    with nuc_insertions.open("w") as f:
        f.write("read_id\tmain\n")
        f.write(insertions)
    # get amino acid insertions from the nextclade.tsv file column "aaInsertions"
    # and write it to aa_inerstions.tsv with the header "read_id" and the{gene_name}
    # TODO: aa_insertions.tsv (in out test data we have no aa insertions.. so hard to test this)
    # for now just make file with header read_id and {gene_name}
    # make a line of each read_id add an empty [] for each gene name
    aa_insertions_f = output_dir / "aa_insertions.tsv"
    # build header
    header = "read_id" + "\t" + "\t".join([f"{gene_name}" for gene_name in gene_names])
    # read from nextclade.tsv file,  colunm is "seqName"
    read_ids = []
    with (input_dir / "nextclade.tsv").open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            read_ids.append(row["seqName"])
    with aa_insertions_f.open("w") as f:
        f.write(header + "\n")
        for read_id in read_ids:
            f.write(read_id + "\t" + "\t".join(["[]" for _ in gene_names]) + "\n")
    # amino acid insertions look like ORF1b:814:D ,S:214:EPE

    with (input_dir / "nextclade.tsv").open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        aa_insertions_data = {
            read_id: {gene: [] for gene in gene_names} for read_id in read_ids
        }
        for row in reader:
            read_id = row["seqName"]
            aa_insertions = row["aaInsertions"]
            if aa_insertions:
                for insertion in aa_insertions.split(","):
                    gene, position, new_amino_acid = insertion.split(":")
                    formatted_insertion = f"{position}:{new_amino_acid}"
                    if gene in aa_insertions_data[read_id]:
                        aa_insertions_data[read_id][gene].append(formatted_insertion)
                    else:
                        logging.warning(
                            f"Gene {gene} not found in gene names for read_id {read_id}"
                        )

    with aa_insertions_f.open("w") as f:
        f.write(header + "\n")
        for read_id in read_ids:
            row_data = [
                ",".join(aa_insertions_data[read_id][gene]) for gene in gene_names
            ]
            f.write(read_id + "\t" + "\t".join(row_data) + "\n")

    # make metadata per read_id
    # get the metadata from the metadata.json file
    # and write it to the metadata.tsv file
    metadata_tsv = output_dir / "metadata.tsv"
    with metadata_file.open() as f:
        metadata = json.load(f)
    # validate that the meatdata keys are the same as defined in the database_config.yaml file as schema /metadata /name
    # if not raise an error
    with database_config.open() as f:
        database_schema = yaml.safe_load(f)
    metadata_keys = set(metadata.keys())
    schema_keys = set([item["name"] for item in database_schema["schema"]["metadata"]])
    # exclude the read_id from the schema keys
    schema_keys.remove("read_id")
    logging.debug(f"Metadata keys: {metadata_keys}")
    logging.debug(f"Schema keys: {schema_keys}")
    if metadata_keys != schema_keys:
        logging.error(
            f"Metadata keys do not match schema keys: {metadata_keys} != {schema_keys}"
        )
        raise ValueError(
            f"Metadata keys do not match schema keys: {metadata_keys} != {schema_keys}"
        )

    with metadata_tsv.open("w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["read_id"] + list(metadata.keys()))
        for read_id in read_ids:
            f.write(
                read_id
                + "\t"
                + "\t".join([str(metadata[key]) for key in metadata.keys()])
                + "\n"
            )

    # get the database_config.yaml file and write it to the nextclade directory
    # copy over the scripts/database_config.yaml file to the nextclade directory
    nextclade_database_config = output_dir / "database_config.yaml"
    with database_config.open() as f:
        database_config = f.read()
    with nextclade_database_config.open("w") as f:
        f.write(database_config)

    # get the reference sequence and write it to the nextclade directory
    # TODO: check if reference genome is really used by silo-input-transformer
    # if so the file should be extracted from the nexclade translation used, see metadata for the reference$
    # for now we just copy the reference genome from the scripts directory
    reference_genome = output_dir / "reference_genomes.json"

    # get the reference genome from the scripts directory
    reference_genome_source = Path("scripts/reference_genomes.json")
    with reference_genome_source.open() as f:
        reference_genome_data = f.read()
    with reference_genome.open("w") as f:
        f.write(reference_genome_data)

    # get unaliged_main.tsv // which is just the same as the nuc_main.fasta file ?
    # copy over the nuc_main.fasta file and name it unaligned_main.tsv
    unaligned_main = output_dir / "unaligned_main.fasta"
    with nuc_main.open() as f:
        nuc_main_data = f.read()
    with unaligned_main.open("w") as f:
        f.write(nuc_main_data)

    logging.info(f"Results saved to: {output_dir}")
    # return the paths to
    #    metadata_fp=metadata_tsv,
    #    database_config_fp=nextclade_database_config,
    #    reference_genomes_fp=reference_genome,

    path_to_files = {
        "metadata_fp": metadata_tsv,
        "database_config_fp": nextclade_database_config,
        "reference_genomes_fp": reference_genome,
    }

    return path_to_files


def transform_to_ndjson(
    sequence_file_directory: Path,
    trafo_config_fp: Path,
    output_dir: Path,
    metadata_fp: Path,
    database_config_fp: Path,
    reference_genomes_fp: Path,
    file_prefixes: dict = {
        "amino_acid_sequence": "gene_",
        "nucleotide_sequence": "nuc_",
        "unaligned_nucleotide_sequence": "unaligned_",
    },
    batch_size: int = 1000,
) -> None:
    """Transforms the sequences to NDJSON format for SILO database."""

    # make a trafo_config.yaml file
    trafo_config = {
        "file_inputs": {
            "metadata": str(metadata_fp),
            "database_config": str(database_config_fp),
            "reference_genomes": str(reference_genomes_fp),
            "sequence_file_directory": str(sequence_file_directory),
        },
        "file_prefixes": file_prefixes,
        "output_dir": str(output_dir),
        "batch_size": batch_size,
    }
    output_dir.mkdir(parents=True, exist_ok=True)  # Ensure the output directory exists
    with trafo_config_fp.open("w") as f:
        yaml.dump(trafo_config, f)
    logging.info(f"Trafo config saved to: {trafo_config_fp}")

    # run the silo_input_transformer with the trafo_config.yaml file
    logging.info(f"Running silo_input_transformer with config: {trafo_config_fp}")
    silo_input_transformer.run_with_config(str(trafo_config_fp))
    logging.info(f"Results saved to: {output_dir}")
    return None


def make_submission_file(result_dir: Path, srLink: str) -> Path:
    """Create a submission file with the given S3 link.

    Args:
        result_dir (Path): The directory to save the submission file.
        srLink (str): The S3 link to include in the submission file.

    Returns:
        Path: The path to the created submission file.
    """
    result_dir_submission = result_dir / "submission"
    result_dir_submission.mkdir(parents=True, exist_ok=True)

    submission_metadata_fp = result_dir_submission / "metadata.tsv"
    with submission_metadata_fp.open("w") as f:
        f.write("submissionId\ts3Link\tversionComment\n")
        f.write(f"001\t{srLink}\t\n")
    logging.info(f"Submission metadata saved to: {submission_metadata_fp}")

    return submission_metadata_fp


def process_directory(
    input_dir: Path,
    sample_id: str,
    batch_id: str,
    result_dir: Path,
    nextclade_reference: str,
    timeline_file: Path,
    primers_file: Path,
    file_name: str = "REF_aln_trim.bam",
    database_config: Path = Path("scripts/database_config.yaml"),
) -> None:
    """Process all files in a given directory.

    Args:
        input_dir (Path): The directory to process. i.e. the directory containing the BAM file.
                          to reach samples/A1_05_2024_10_08/20241024_2411515907/alignments/
        result_dir (Path): The directory to save the results.
        nextclade_reference (str): The reference to use for nextclade.
        timeline_file (Path): The timeline file to cross-reference the metadata.
        primers_file (Path): The primers file to cross-reference the metadata.
        file_name (str): The name of the file to process

    Returns:
        None (writes results to the result_dir)
    """

    # TODO: absolb all these intermediary files into a temporary directory

    # check that one was given a directory and not a file and it exists
    if not input_dir.is_dir():
        logging.error(f"Input directory not found, is it a directory?: {input_dir}")
        raise FileNotFoundError(f"Directory not found: {input_dir}")

    logging.info(f"Processing directory: {input_dir}")
    logging.info(f"Assuming the input file is: {file_name}")
    # check that the file exists and also it's .bai file
    sample_fp = input_dir / file_name
    if not sample_fp.exists():
        logging.error(f"Input file not found: {sample_fp}")
        raise FileNotFoundError(f"Input file not found: {sample_fp}")

    ##### Get Sample and Batch metadata and write to a file #####
    metadata = get_metadata(sample_id, batch_id, timeline_file, primers_file)
    # add nextclade reference to metadata
    metadata["nextclade_reference"] = nextclade_reference
    metadata_file = result_dir / "metadata.json"
    result_dir.mkdir(parents=True, exist_ok=True)
    with metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)
    logging.info(f"Metadata saved to: {metadata_file}")

    ##### Convert BAM to SAM #####
    logging.info(f"Converting BAM to SAM")
    bam_file = sample_fp
    sam_data = bam_to_sam(bam_file)

    ##### Process SAM to FASTA #####
    logging.info(f"Processing SAM to FASTA (pair, merge, and normalize reads)")
    fasta_file = result_dir / "reads.fasta"
    insertions_file = result_dir / "insertions.txt"
    pair_normalize_reads(sam_data, fasta_file, insertions_file)

    ##### Translate nucleotides to amino acids #####
    logging.info(f"Aliging and translating sequences")
    results_dir_translated = result_dir / "translated"
    translate([fasta_file], results_dir_translated, nextclade_reference)

    logging.info(f"Results saved to: {results_dir_translated}")

    ##### Wrangle to Nextclade format // silo_input_transformer inputs #####
    result_dir_wrangled = result_dir / "wrangled"
    path_to_files = wrangle_for_transformer(
        input_dir=results_dir_translated,
        output_dir=result_dir_wrangled,
        fasta_file=fasta_file,
        insertions_file=insertions_file,
        metadata_file=metadata_file,
        database_config=database_config,
    )

    ###### Transform to NDJSON ######
    result_dir_transformed = result_dir / "transformed"
    logging.debug(f"Transforming to NDJSON")
    logging.debug(f"sequence_file_directory: {result_dir_wrangled}")
    logging.debug(f"output_dir: {result_dir_transformed}")
    transform_to_ndjson(
        sequence_file_directory=result_dir_wrangled,
        trafo_config_fp=result_dir / "trafo_config.yaml",
        output_dir=result_dir_transformed,
        metadata_fp=path_to_files["metadata_fp"],
        database_config_fp=path_to_files["database_config_fp"],
        reference_genomes_fp=path_to_files["reference_genomes_fp"],
    )

    #####   Compress & Upload to S3  #####
    file_to_upload = result_dir_transformed / "silo_input.ndjson"
    compressed_file = result_dir_transformed / "silo_input.ndjson.bz2"
    logging.info(f"Compressing file: {file_to_upload}")
    compress_bz2(file_to_upload, compressed_file)

    #  Upload as generate a file name for the submission file, i.e. use the SAMPLE_ID
    logging.info(f"Uploading to S3: {compressed_file}")
    s3_file_name = f"{sample_id}.ndjson.bz2"
    s3_bucket = "sr2silo01"
    s3_link = f"s3://{s3_bucket}/{s3_file_name}"
    upload_file_to_s3(compressed_file, s3_bucket, s3_file_name)

    ##### Submit S3 reference to SILO #####
    logging.info(f"Submitting to Loculus")
    input_fp = make_submission_file(result_dir, s3_link)
    username = "testuser"
    password = "testuser"
    group_id = 1
    submit(input_fp, username, password, group_id)


@click.command()
@click.option("--sample_dir", envvar="SAMPLE_DIR", help="Path to the sample directory.")
@click.option("--sample_id", envvar="SAMPLE_ID", help="sample_id to use for metadata.")
@click.option("--batch_id", envvar="BATCH_ID", help="batch_id to use for metadata.")
@click.option(
    "--result_dir", envvar="RESULTS_DIR", help="Path to the results directory."
)
@click.option(
    "--timeline_file", envvar="TIMELINE_FILE", help="Path to the timeline file."
)
@click.option("--primer_file", envvar="PRIMER_FILE", help="Path to the primers file.")
@click.option(
    "--nextclade_reference",
    envvar="NEXTCLADE_REFERENCE",
    default="sars-cov-2",
    help="Nextclade reference.",
)
def main(
    sample_dir,
    sample_id,
    batch_id,
    result_dir,
    timeline_file,
    primer_file,
    nextclade_reference,
    database_config: Path = Path("scripts/database_config.yaml"),
):
    """Process a sample directory."""
    logging.info(f"Processing sample directory: {sample_dir}")
    logging.info(f"Saving results to: {result_dir}")
    logging.info(f"Using timeline file: {timeline_file}")
    logging.info(f"Using primers file: {primer_file}")
    logging.info(f"Using Nextclade reference: {nextclade_reference}")
    logging.info(f"Using sample_id: {sample_id}")
    logging.info(f"Using batch_id: {batch_id}")
    logging.info(f"Using database_config: {database_config}")

    process_directory(
        input_dir=Path("sample"),
        sample_id=sample_id,
        batch_id=batch_id,
        result_dir=Path("results"),
        timeline_file=Path("timeline.tsv"),
        primers_file=Path("primers.yaml"),
        nextclade_reference=nextclade_reference,
        database_config=Path("scripts/database_config.yaml"),
    )


if __name__ == "__main__":
    main()
