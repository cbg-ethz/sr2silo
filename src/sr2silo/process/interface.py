"""Implements Genomic Data Structures"""

from __future__ import annotations

import json
import logging
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel, ValidationError

from sr2silo.silo_aligned_read import AlignedReadSchema, ReadMetadata

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


class NucInsertion:
    """A nuclotide insertion."""

    def __init__(self, position: int, sequence: str):
        """Initialize with a position and a sequence."""
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        """toString method."""
        return f"{self.position} : {self.sequence}"


class AAInsertion:
    """An amino acid insertion."""

    def __init__(self, position: int, sequence: str):
        """Initialize with a position and a sequence."""
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        """toString method."""
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
        metadata: Optional[Union[Dict[str, str], ReadMetadata]] = None,
    ):
        """Initialize with a AlignedRead object."""
        self.read_id = read_id
        self.unaligned_nucleotide_sequences = unaligned_nucleotide_sequences
        self.aligned_nucleotide_sequences = aligned_nucleotide_sequences
        self.nucleotide_insertions = nucleotide_insertions
        self.amino_acid_insertions = amino_acid_insertions
        self.aligned_amino_acid_sequences = aligned_amino_acid_sequences
        self.metadata = metadata
        self._validate_types()

    def _validate_types(self):
        """Validate the types of the attributes."""
        if not isinstance(self.read_id, str):
            raise TypeError(f"read_id must be a str, got {type(self.read_id).__name__}")
        if not isinstance(self.unaligned_nucleotide_sequences, str):
            raise TypeError(
                f"unaligned_nucleotide_sequences must be a str, got "
                f"{type(self.unaligned_nucleotide_sequences).__name__}"
            )
        if not isinstance(self.aligned_nucleotide_sequences, str):
            raise TypeError(
                f"aligned_nucleotide_sequences must be a str, got "
                f"{type(self.aligned_nucleotide_sequences).__name__}"
            )
        if not all(isinstance(i, NucInsertion) for i in self.nucleotide_insertions):
            raise TypeError(
                "All items in nucleotide_insertions must be NucInsertion instances"
            )
        if not isinstance(self.amino_acid_insertions, AAInsertionSet):
            raise TypeError(
                f"amino_acid_insertions must be an AAInsertionSet, got "
                f"{type(self.amino_acid_insertions).__name__}"
            )
        if not isinstance(self.aligned_amino_acid_sequences, AASequenceSet):
            raise TypeError(
                f"aligned_amino_acid_sequences must be a dict, got "
                f"{type(self.aligned_amino_acid_sequences).__name__}"
            )
        # TODO: rework what metadata is allowed to be
        if self.metadata is not None and not isinstance(
            self.metadata, (dict, ReadMetadata)
        ):
            raise TypeError(
                "metadata must be a dict or ReadMetadata, "
                f"got {type(self.metadata).__name__}"
            )

    def set_nuc_insertion(self, nuc_insertion: NucInsertion):
        """Append a nucleotide insertion to the list of nucleotide insertions."""
        self.nucleotide_insertions.append(nuc_insertion)

    def get_amino_acid_insertions(self) -> AAInsertionSet:
        """Return the amino acid insertions."""
        return self.amino_acid_insertions

    def set_metadata(self, metadata: Union[Dict[str, str], BaseModel]):
        """Set the metadata. If a BaseModel is provided, convert it to dict."""
        if isinstance(metadata, BaseModel):
            self.metadata = metadata.model_dump()
        else:
            self.metadata = metadata

    def get_metadata(self) -> Optional[Union[Dict[str, str], BaseModel]]:
        """Return the metadata."""
        return self.metadata

    def to_dict(self) -> Dict[str, Any]:
        """Return a dictionary / json representation of the object."""
        formatted_nuc_ins = [
            f"{ins.position} : {ins.sequence}" for ins in self.nucleotide_insertions
        ]
        json_representation = {
            "readId": self.read_id,
            "nucleotideInsertions": {
                "main": formatted_nuc_ins,
            },
            "aminoAcidInsertions": self.amino_acid_insertions.to_dict(),
            "alignedNucleotideSequences": {
                "main": self.aligned_nucleotide_sequences,
            },
            "unalignedNucleotideSequences": {
                "main": self.unaligned_nucleotide_sequences,
            },
            "alignedAminoAcidSequences": self.aligned_amino_acid_sequences.to_dict(),
        }
        if self.metadata:
            if isinstance(self.metadata, dict):
                self.metadata["read_id"] = self.read_id
            json_representation["metadata"] = self.metadata
        return json_representation

    def to_silo_json(self) -> str:
        """
        Validate the aligned read dict using a pydantic schema and print a
        nicely formatted JSON string conforming to the DB requirements.
        """
        try:
            schema = AlignedReadSchema(**self.to_dict())
            return schema.model_dump_json(indent=2, exclude_none=True)
        except ValidationError as e:
            raise e

    def __str__(self) -> str:
        """toString method as pretty JSON string."""
        return json.dumps(self.to_dict(), indent=2)

    @staticmethod
    def from_str(data: str) -> AlignedRead:
        """Create an AlignedRead object from a string."""
        data = data.strip()  # Remove extra whitespace

        # Parse the json data to a dict
        json_data = json.loads(data)
        json_data = json.loads(json_data)

        read_id = json_data["readId"]
        unaligned_nucleotide_sequences = json_data["unalignedNucleotideSequences"][
            "main"
        ]
        aligned_nucleotide_sequences = json_data["alignedNucleotideSequences"]["main"]
        nucleotide_insertions = []
        if json_data["nucleotideInsertions"]["main"]:
            nucleotide_insertions = [
                NucInsertion(int(ins.split(" : ")[0]), ins.split(" : ")[1])
                for ins in json_data["nucleotideInsertions"]["main"]
            ]
        amino_acid_insertions = AAInsertionSet.from_dict(
            json_data["aminoAcidInsertions"]
        )
        aligned_amino_acid_sequences = AASequenceSet.from_dict(
            json_data["alignedAminoAcidSequences"]
        )

        # validate all the arguments are of the correct type
        assert isinstance(read_id, str)
        assert isinstance(unaligned_nucleotide_sequences, str)
        assert isinstance(aligned_nucleotide_sequences, str)
        assert all(isinstance(i, NucInsertion) for i in nucleotide_insertions)
        assert isinstance(amino_acid_insertions, AAInsertionSet)
        assert isinstance(aligned_amino_acid_sequences, AASequenceSet)

        try:
            return AlignedRead(
                read_id,
                unaligned_nucleotide_sequences,
                aligned_nucleotide_sequences,
                nucleotide_insertions,
                amino_acid_insertions,
                aligned_amino_acid_sequences,
            )
        except TypeError as e:
            logging.error(
                "Error constructing AlignedRead with data: " + repr(json_data)
            )
            raise e


class GeneName:
    """Class to represent a gene name in its short form."""

    def __init__(self, name: str):
        """Initialize with a gene name."""
        self.name = name

    def __str__(self) -> str:
        """toString method."""
        return self.name


class Gene:
    def __init__(self, gene_name: GeneName, gene_length: int):
        """Initialize with a gene name and a gene length."""
        self.name = gene_name
        self.gene_length = gene_length

    def to_dict(self) -> Dict[str, GeneName | int]:
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
        self.genes[str(gene_name)].gene_length = gene_length

    def get_gene(self, gene_name: GeneName) -> Gene:
        """Return a gene by name."""
        return self.genes[str(gene_name)]

    def get_gene_length(self, gene_name: GeneName) -> int:
        """Return the length of a gene."""
        return self.genes[str(gene_name)].gene_length

    def get_gene_name_list(self) -> List[GeneName]:
        """Return a list of genes."""
        return list(self.genes.keys())

    def to_dict(self) -> Dict[str, Dict[str, Any]]:
        """Return a dictionary with gene names as keys and gene
        length as values."""
        return {
            str(k): {"gene_name": str(v.name), "gene_length": v.gene_length}
            for k, v in self.genes.items()
        }

    def __str__(self) -> str:
        return str(self.to_dict())


class AAInsertionSet:
    """Class to represent the set of amino acid insertions for a full set of genes."""

    def __init__(self, genes: List[GeneName]):
        """Initialize with an empty set of insertions for each gene."""
        self.aa_insertions = {gene: [] for gene in genes}

    def set_insertions_for_gene(
        self, gene_name: GeneName, aa_insertions: List[AAInsertion]
    ):
        """Set the amino acid insertions for a particular gene."""
        self.aa_insertions[gene_name] = aa_insertions

    def to_dict(self) -> dict:
        """Return a dictionary with gene names as keys."""
        return {
            str(gene): [f"{ins.position} : {ins.sequence}" for ins in ins_per_gene]
            for gene, ins_per_gene in self.aa_insertions.items()
        }

    def __str__(self) -> str:
        return str(self.to_dict())

    @staticmethod
    def from_dict(data: Dict) -> AAInsertionSet:
        """Create an AAInsertionSet object from a dictionary."""
        aa_insertions = AAInsertionSet([])
        for gene_name, ins_list in data.items():
            aa_insertions.aa_insertions[gene_name] = [
                AAInsertion(int(ins.split(" : ")[0]), ins.split(" : ")[1])
                for ins in ins_list
            ]
        return aa_insertions


class AASequenceSet:
    """Class to represent the set of amino acid sequences for a full set of genes."""

    def __init__(self, genes: List[GeneName]):
        """Initialize with an empty sequence for each gene"""
        self.sequences = {gene: "" for gene in genes}
        self.genes = genes

    def set_sequence(self, gene_name: GeneName, aa_sequence: str):
        """Set the amino acid sequence for a particular gene."""
        self.sequences[gene_name] = aa_sequence

    def to_dict(self) -> dict:
        """Return a dictionary with gene names as keys"""
        return {str(gene): seq for gene, seq in self.sequences.items()}

    @staticmethod
    def from_dict(data: dict) -> AASequenceSet:
        """Create an AASequenceSet object from a dictionary."""
        aa_sequences = AASequenceSet([])
        for gene_name, seq in data.items():
            aa_sequences.sequences[gene_name] = seq
        return aa_sequences

    def __str__(self) -> str:
        return str(self.to_dict())
