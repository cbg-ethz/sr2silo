"""Implements Genomic Data Structures"""


from __future__ import annotations

import json

class NucInsertion:
    """A nuclotide insertion."""
    def __init__(self, position: int, sequence: str):
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        return f"{self.position} : {self.sequence}"


class AAInsertion:
    """An amino acid insertion."""
    def __init__(self, position: int, sequence: str):
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        return f"{self.position} : {self.sequence}"


class AlignedRead:
    """Class to represent an aligned read."""
    def __init__(
        self,
        read_id: str,
        unaligned_nucleotide_sequences: str,
        aligned_nucleotide_sequences: str,
        nucleotide_insertions: List[NucInsertion],
        amino_acid_insertions: AAInsertionSet,
        aligned_amino_acid_sequences: AASequenceSet,
    ):
        self.read_id = read_id
        self.unaligned_nucleotide_sequences = unaligned_nucleotide_sequences
        self.aligned_nucleotide_sequences = aligned_nucleotide_sequences
        self.nucleotide_insertions = nucleotide_insertions
        self.amino_acid_insertions = amino_acid_insertions
        self.aligned_amino_acid_sequences = aligned_amino_acid_sequences

        self._validate_types()

    def _validate_types(self):
        if not isinstance(self.read_id, str):
            raise TypeError(f"read_id must be a str, got {type(self.read_id).__name__}")
        if not isinstance(self.unaligned_nucleotide_sequences, str):
            raise TypeError(f"unaligned_nucleotide_sequences must be a str, got {type(self.unaligned_nucleotide_sequences).__name__}")
        if not isinstance(self.aligned_nucleotide_sequences, str):
            raise TypeError(f"aligned_nucleotide_sequences must be a str, got {type(self.aligned_nucleotide_sequences).__name__}")
        if not all(isinstance(i, NucInsertion) for i in self.nucleotide_insertions):
            raise TypeError("All items in nucleotide_insertions must be NucInsertion instances")
        if not isinstance(self.amino_acid_insertions, AAInsertionSet):
            raise TypeError(f"amino_acid_insertions must be an AAInsertionSet, got {type(self.amino_acid_insertions).__name__}")
        if not isinstance(self.aligned_amino_acid_sequences, AASequenceSet):
            raise TypeError(f"aligned_amino_acid_sequences must be a dict, got {type(self.aligned_amino_acid_sequences).__name__}")

    def set_nuc_insertion(self, nuc_insertion: NucInsertion):
        """Append a nucleotide insertion to the list of nucleotide insertions."""
        self.nucleotide_insertions.append(nuc_insertion)

    def to_dict(self) -> Dict[str, any]:  # noqa
        return {
            "read_id": self.read_id,
            "unaligned_nucleotide_sequences": self.unaligned_nucleotide_sequences,
            "aligned_nucleotide_sequences": self.aligned_nucleotide_sequences,
            "nucleotide_insertions": self.nucleotide_insertions,
            "amino_acid_insertions": self.amino_acid_insertions,
            "aligned_amino_acid_sequences": self.aligned_amino_acid_sequences,
        }

    def get_amino_acid_insertions(self) -> AAInsertionSet:
        return self.amino_acid_insertions

    def __str__(self) -> str:
        json_representation = {
            "readId": self.read_id,
            "nucleotideInsertions": {
                "main": [str(ins) for ins in self.nucleotide_insertions],
            },
            "aminoAcidInsertions": self.amino_acid_insertions.__str__(),
            "alignedNucleotideSequences": {
                "main": self.aligned_nucleotide_sequences,
            },
            "unalignedNucleotideSequences": {
                 "main": self.unaligned_nucleotide_sequences,
            },
             "alignedAminoAcidSequences": self.aligned_amino_acid_sequences.__str__(),
        }

        return str(json_representation)


class GeneName:
    """Class to represent a gene name in its short form."""
    def __init__(self, name: str):
        self.name = name

    def __str__(self) -> str:
        return self.name


class Gene:
    def __init__(self, gene_name: GeneName, gene_length: int):
        self.name = gene_name
        self.gene_length = gene_length

    def to_dict(self) -> Dict[str, int | str]:
        return {
            "gene_name": self.name,
            "gene_length": self.gene_length,
        }


class GeneSet:
    """Class to represent a set of genes for a pathogen"""
    def __init__(self, genes: List[Gene]):
        """Initialize with a list of genes."""
        self.genes = {}
        for gene in genes:
            self.genes[str(gene.name)] = gene

    def set_gene_length(self, gene_name: GeneName, gene_length: int):
        """Set the length of a gene."""
        self.genes[str(name)].gene_length = gene_length

    def get_gene(self, gene_name: GeneName) -> Gene:
        """Return a gene by name."""
        return self.genes[str(gene_name)]

    def get_gene_length(self, gene_name: GeneName) -> int:
        """Return the length of a gene."""
        return self.genes[str(gene_name)].gene_length

    def get_gene_name_list(self) -> List[GeneName]:
        """Return a list of genes."""
        return list(self.genes.keys())

    def to_dict(self) -> Dict[str, Dict[str, int | str]]:
        """Return a dictionary with gene names as keys and gene
           leght as values."""
        return ({str(k): v} for k,v in self.genes.items())

    def __str__(self) -> str:
        return str(self.to_dict())

class AAInsertionSet:
    """Class to represent the set of amino acid insertions for a full set of genes."""

    def __init__(self, genes: List[GeneName]):
        """Initialize with an empty set of insertions for each gene."""
        self.aa_insertions = {str(gene): [] for gene in genes}

    def set_insertions_for_gene(self, gene_name: GeneName, aa_insertions: List[AAInsertion]):
        """Set the amino acid insertions for a particular gene."""
        self.aa_insertions[gene_name] = aa_insertions

    def to_dict(self) -> dict:
        """Return a dictionary with gene names as keys."""
        return {str(gene): [str(ins) for ins in ins_per_gene] for gene, ins_per_gene in self.aa_insertions.items() if isinstance(ins_per_gene, list)}

    def __str__(self) -> str:
        return str(self.to_dict())


class AASequenceSet:
    """Class to represent the set of amino acid sequences for a full set of genes."""

    def __init__(self, genes: List[GeneName]):
        """ Initialize with an empty sequence for each gene"""
        self.sequences = {gene: "" for gene in genes}
        self.genes = genes

    def set_sequence(self, gene_name: str, aa_sequence: str):
        """Set the amino acid sequence for a particular gene."""
        self.sequences[gene_name] = aa_sequence

    def to_dict(self) -> dict:
        """Return a dictionary with gene names as keys"""
        return {str(gene): seq for gene, seq in self.sequences.items()}

    def __str__(self) -> str:
        return str(self.to_dict())

