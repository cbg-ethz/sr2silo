"""Wrangles the data already processed from V-Pipe and trandlated to amino acids
   to a format that SILO and LAPIS accept.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import yaml


# TODO: to be refactored once the silo_import_transformer is moved to SILO
def wrangle_for_transformer(
    input_dir: Path,
    output_dir: Path,
    fasta_file: Path,
    insertions_file: Path,
    metadata_file: Path,
    database_config_fp: Path,
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
        metadata_file (Path): The metadata json file containing the per
                              sequencing run metadata, that is copied for each
                              read_id in the output.
        database_config_fp (Path): The database configuration file containing the
                                schema for the metadata, has to match the
                                metadata keys in the metadata.json file but the
                                read_id key is not in the schema.

    Returns:
        dict[str, Path]: A dictionary containing the paths to the files
                         created during the wrangling process.

                         metadata_fp: metadata.tsv
                         database_config_new_fp: database_config.yaml
                         reference_genomes_fp: reference_genomes.json
    """

    logging.info("Wrangling sequences to Nextclade format")
    # make a directory for the nextclade-like results
    output_dir.mkdir(parents=True, exist_ok=True)
    # copy over the nucleotide sequence file and name nuc_main.fasta
    nuc_main = output_dir / "nuc_main.fasta"
    with fasta_file.open() as f:
        nuc_main_data = f.read()
    with nuc_main.open("w") as f:
        f.write(nuc_main_data)
    # copy over the amino acids sequences and name them gene_ and keep the file
    #  ending and last part of the name
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
    # copy over the nucoletide insertions file and name it
    # nucolotide_insertions.tsv
    # add a header of "read_id  | main" to the file
    nuc_insertions = output_dir / "nucleotide_insertions.tsv"
    with insertions_file.open() as f:
        insertions = f.read()
    with nuc_insertions.open("w") as f:
        f.write("read_id\tmain\n")
        f.write(insertions)
    # get amino acid insertions from the nextclade.tsv file column "aaInsertions"
    # and write it to aa_inerstions.tsv with the header "read_id" and the{gene_name}
    # TODO: aa_insertions.tsv (in out test data we have no aa insertions..
    # so hard to test this)
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
    # validate that the meatdata keys are the same as defined in the
    # database_config.yaml file as schema /metadata /name
    # if not raise an error
    with database_config_fp.open() as f:
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
    with database_config_fp.open() as f:
        database_config = f.read()
    with nextclade_database_config.open("w") as f:
        f.write(database_config)

    # get the reference sequence and write it to the nextclade directory
    # TODO: check if reference genome is really used by silo-input-transformer
    # if so the file should be extracted from the nexclade translation used,
    # see metadata for the reference$
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

    path_to_files = {
        "metadata_fp": metadata_tsv,
        "database_config_new_fp": nextclade_database_config,
        "reference_genomes_fp": reference_genome,
    }

    return path_to_files
