"""This module contains the logic to translate consensus nuclotides to
    amino acid sequences."""

from __future__ import annotations

import subprocess


def translate(input_file: str, result_dir: str, dataset_dir: str):
    """Translate consensus nucleotides to amino acid sequences."""

    # first get the test dataset from the gff3 file
    subprocess.run(
        [
            "nextclade",
            "dataset",
            "get",
            "--name",
            "nextstrain/sars-cov-2/wuhan-hu-1/orfs",
            "--output-dir",
            "data/sars-cov-2",
        ],
        check=True,
    )

    # then replace the sequences.fasta in the
    # data/sars-cov-2 with the sequences.fasta from the input file
    subprocess.run(["cp", input_file, "data/sars-cov-2/sequences.fasta"], check=True)

    # then run the nextclade run command
    subprocess.run(
        [
            "nextclade",
            "run",
            "--input-dataset",
            "data/sars-cov-2",
            "--output-all=output/",
            "data/sars-cov-2/sequences.fasta",
        ],
        check=True,
    )

    # # This is how the loculus guys do its
    # command = [
    #     "nextclade3",
    #     "run",
    #     f"--output-all={result_dir_seg}",
    #     f"--input-dataset={dataset_dir_seg}",
    #     f"--output-translations={result_dir_seg}/nextclade.cds_translation.{{cds}}.fasta",
    #     "--jobs=1",
    #     "--",
    #     input_file,
    # ]


#     nextclade run \
#    --verbose \
#    --include-reference \
#    --in-order \
#    --input-dataset=data/sars-cov-2 \
#    --input-ref=data/sars-cov-2/reference.fasta \
#    --input-annotation=data/sars-cov-2/genome_annotation.gff3 \
#    --cds-selection=E,M,N,ORF1a,ORF1b,ORF3a,ORF6,ORF7a,ORF7b,ORF8,ORF9b,S \
#    --input-tree=data/sars-cov-2/tree.json \
#    --input-pathogen-json=data/sars-cov-2/pathogen.json \
#    --output-fasta=output/nextclade.aligned.fasta.gz \
#    --output-json=output/nextclade.json \
#    --output-ndjson=output/nextclade.ndjson \
#    --output-csv=output/nextclade.csv \
#    --output-tsv=output/nextclade.tsv \
#    --output-tree=output/nextclade.auspice.json \
#    --output-tree-nwk=output/nextclade.tree.nwk \
#    --output-translations=output/nextclade_CDS_{cds}.translation.fasta.zst \
#    data/sars-cov-2/sequences.fasta \
#    my_sequences1.fasta.gz \
#    my_sequences2.fasta.xz

#     subprocess.run(command, check=True)


#     nextstrain/sars-cov-2/XBB
