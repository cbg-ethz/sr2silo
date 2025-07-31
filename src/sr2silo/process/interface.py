"""Implements Genomic Data Structures"""

from __future__ import annotations

import json
import logging
from typing import Any, Dict, List, Optional, Union

from pydantic import BaseModel

from sr2silo.silo_read_schema import AlignedReadSchema, ReadMetadata

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


class Insertion:
    """Base class for insertions."""

    def __init__(self, position: int, sequence: str):
        """Initialize with a position and a sequence."""
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        """toString method."""
        return f"{self.position}:{self.sequence}"

    def __eq__(self, other) -> bool:
        """Compare two Insertion objects by their values.
        This compares position and sequence values,
        ensuring a value-based comparison rather than identity-based.

        Also ensures that the objects are of the exact same class,
        not just derived from the same base class.
        """
        if not isinstance(other, Insertion):
            return False
        # Check if the objects are of the same specific class (not just Insertion)
        if not isinstance(other, self.__class__) or not isinstance(
            self, other.__class__
        ):
            return False
        return self.position == other.position and self.sequence == other.sequence


class NucInsertion(Insertion):
    """A nuclotide insertion."""

    pass


class AAInsertion(Insertion):
    """An amino acid insertion."""

    pass


class AlignedRead:
    """Class to represent an aligned read."""

    __slots__ = [
        "read_id",
        "unaligned_nucleotide_sequence",
        "aligned_nucleotide_sequence",
        "aligned_nucleotide_sequence_offset",
        "nucleotide_insertions",
        "amino_acid_insertions",
        "aligned_amino_acid_sequences",
        "metadata",
    ]

    def __init__(
        self,
        read_id: str,
        unaligned_nucleotide_sequence: str,
        aligned_nucleotide_sequence: str,
        aligned_nucleotide_sequence_offset: int,
        nucleotide_insertions: List[NucInsertion],
        amino_acid_insertions: AAInsertionSet,
        aligned_amino_acid_sequences: AASequenceSet,
        metadata: Optional[Union[Dict[str, str], ReadMetadata]] = None,
    ):
        """Initialize with a AlignedRead object."""
        self.read_id = read_id
        self.unaligned_nucleotide_sequence = unaligned_nucleotide_sequence
        self.aligned_nucleotide_sequence = aligned_nucleotide_sequence
        self.nucleotide_insertions = nucleotide_insertions
        self.amino_acid_insertions = amino_acid_insertions
        self.aligned_amino_acid_sequences = aligned_amino_acid_sequences
        self.aligned_nucleotide_sequence_offset = aligned_nucleotide_sequence_offset
        if metadata:
            self.set_metadata(metadata)
        else:
            self.metadata = None
        self._validate_types()

    def _validate_types(self):
        """Validate the types of the attributes."""
        if not isinstance(self.read_id, str):
            raise TypeError(f"read_id must be a str, got {type(self.read_id).__name__}")
        if not isinstance(self.unaligned_nucleotide_sequence, str):
            raise TypeError(
                f"unaligned_nucleotide_sequence must be a str, got "
                f"{type(self.unaligned_nucleotide_sequence).__name__}"
            )
        if not isinstance(self.aligned_nucleotide_sequence, str):
            raise TypeError(
                f"aligned_nucleotide_sequence must be a str, got "
                f"{type(self.aligned_nucleotide_sequence).__name__}"
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
        if isinstance(metadata, ReadMetadata):
            self.metadata = metadata.model_dump()
        else:
            self.metadata = metadata

    def get_metadata(self) -> Optional[Union[Dict[str, str], ReadMetadata]]:
        """Return the metadata."""
        if self.metadata:
            if isinstance(self.metadata, dict):
                return self.metadata
            else:
                return self.metadata.dict()
        return None

    def to_dict(self) -> Dict[str, Any]:
        """Return a dictionary / json representation of the object."""
        formatted_nuc_ins = [
            f"{ins.position}:{ins.sequence}" for ins in self.nucleotide_insertions
        ]
        aln_gene_list = aa_sequence_set_and_insertions_to_aligned_genes(
            self.aligned_amino_acid_sequences, self.amino_acid_insertions
        )

        json_representation = {}

        # Add metadata at the start if available
        if self.metadata:
            metadata_dict = (
                self.metadata.model_dump()
                if isinstance(self.metadata, BaseModel)
                else self.metadata
            )
            for key, value in metadata_dict.items():
                if isinstance(value, str):
                    json_representation[key] = value
                elif isinstance(value, (int, float)):
                    json_representation[key] = str(value)
                else:
                    logging.warning(
                        f"Metadata value for {key} is not a string or number: {value}"
                    )

        # Add main section
        json_representation["main"] = {
            "insertions": formatted_nuc_ins,
            "sequence": self.aligned_nucleotide_sequence,
            "offset": self.aligned_nucleotide_sequence_offset,
        }

        # Add each gene to the JSON representation
        for aligned_gene in aln_gene_list:
            # Check if gene has no sequence - if so, represent as null
            # Note: Insertions without a context sequence are non-sense,
            # so we only need to check if the sequence is empty
            if not aligned_gene.sequence:
                json_representation[aligned_gene.gene_name.name] = None
            else:
                formatted_aa_ins = [str(ins) for ins in aligned_gene.insertions]
                json_representation[aligned_gene.gene_name.name] = {
                    "sequence": aligned_gene.sequence,
                    "offset": aligned_gene.offset,
                    "insertions": formatted_aa_ins,
                }

        return json_representation

    def to_silo_json(self, indent: bool = False) -> str:
        """
        Validate the aligned read dict using the new SILO schema and return a
        nicely formatted JSON string conforming to the DB requirements.

        Args:
            indent: Whether to indent the JSON string, True for pretty print.
        """
        # Get the dictionary representation
        data_dict = self.to_dict()

        try:
            # Validate using the new AlignedReadSchema
            schema = AlignedReadSchema(**data_dict)

            return schema.model_dump_json(
                indent=2 if indent else None, exclude_none=False
            )
        except Exception as e:
            # If validation fails, log the error and return the raw JSON
            logging.error(f"Schema validation failed: {e}")
            logging.error(f"Data that failed validation: {data_dict}")

            # Return unvalidated JSON for debugging
            return json.dumps(
                data_dict, indent=2 if indent else None, ensure_ascii=False
            )

    def __str__(self) -> str:
        """toString method as pretty JSON string."""
        return json.dumps(self.to_dict(), indent=2)


class GeneName:
    """Class to represent a gene name in its short form."""

    def __init__(self, name: str):
        """Initialize with a gene name."""
        self.name = name

    def __str__(self) -> str:
        """toString method."""
        return self.name


class Gene:
    """Class to represent a gene with a name and a length."""

    def __init__(self, gene_name: GeneName, gene_length: int):
        """Initialize with a gene name and a gene length."""
        self.name = gene_name
        self.gene_length = gene_length

    def to_dict(self) -> Dict[str, GeneName | int]:
        """Return a dictionary with gene name and gene length."""
        return {
            "gene_name": self.name,
            "gene_length": self.gene_length,
        }

    def __str__(self) -> str:
        gene_str = f"gene_name: {self.name.name}, gene_length: {self.gene_length}"
        return "{" + gene_str + "}"


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
            str(gene): [f"{ins.position}:{ins.sequence}" for ins in ins_per_gene]
            for gene, ins_per_gene in self.aa_insertions.items()
        }

    def __str__(self) -> str:
        return str(self.to_dict())

    def __eq__(self, other) -> bool:
        """Compare two AAInsertionSet objects by their values.

        This compares the dictionaries generated by to_dict() method,
        ensuring a value-based comparison rather than identity-based.
        """
        if not isinstance(other, AAInsertionSet):
            return False

        # Convert both objects to dictionaries for comparison
        self_dict = self.to_dict()
        other_dict = other.to_dict()

        # Check if dictionaries have the same keys
        if set(self_dict.keys()) != set(other_dict.keys()):
            return False

        # For each key, check if the lists have the same elements in the same order
        for key in self_dict:
            if sorted(self_dict[key]) != sorted(other_dict[key]):
                return False

        return True

    @staticmethod
    def from_dict(data: Dict) -> AAInsertionSet:
        """Create an AAInsertionSet object from a dictionary."""
        aa_insertions = AAInsertionSet([])
        for gene_name, ins_list in data.items():
            aa_insertions.aa_insertions[gene_name] = [
                AAInsertion(int(ins.split(":")[0]), ins.split(":")[1])
                for ins in ins_list
            ]
        return aa_insertions


class AASequenceSet:
    """Class to represent the set of amino acid sequences for a full set of genes."""

    def __init__(self, genes: List[GeneName]):
        """Initialize with an empty sequence for each gene"""
        self.sequences = {gene: "" for gene in genes}
        self.offsets = {gene: 0 for gene in genes}
        self.genes = genes

    def set_sequence(self, gene_name: GeneName, aa_sequence: str, offset: int = 0):
        """Set the amino acid sequence for a particular gene."""
        self.sequences[gene_name] = aa_sequence
        self.offsets[gene_name] = offset

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


class AlignedGene:
    """Class to represent an aligned gene with its sequence, insertions, and offset."""

    def __init__(
        self,
        gene_name: GeneName,
        sequence: str,
        offset: int,
        insertions: List[AAInsertion],
    ):
        """Initialize with gene name, sequence, offset, and insertions."""
        self.gene_name = gene_name
        self.sequence = sequence
        self.offset = offset
        self.insertions = insertions

    def __str__(self) -> str:
        """String representation of the aligned gene."""
        return (
            f"AlignedGene(gene_name={self.gene_name}, "
            f"sequence_length={len(self.sequence)}, offset={self.offset}, "
            f"insertions_count={len(self.insertions)})"
        )

    def __eq__(self, other) -> bool:
        """Compare two AlignedGene objects by their values."""
        if not isinstance(other, AlignedGene):
            return False
        return (
            self.gene_name.name == other.gene_name.name
            and self.sequence == other.sequence
            and self.offset == other.offset
            and self.insertions == other.insertions
        )


def aa_sequence_set_and_insertions_to_aligned_genes(
    aa_sequence_set: AASequenceSet,
    aa_insertion_set: AAInsertionSet,
) -> List[AlignedGene]:
    """
    Convert an AASequenceSet and AAInsertionSet to a list of AlignedGene objects.

    Args:
        aa_sequence_set: The amino acid sequences for each gene
        aa_insertion_set: The amino acid insertions for each gene

    Returns:
        List of AlignedGene objects
    """
    aligned_genes = []

    # Get all gene names from the sequence set
    for gene_name_str, sequence in aa_sequence_set.sequences.items():
        # Create GeneName object if gene_name_str is a string
        if isinstance(gene_name_str, str):
            gene_name = GeneName(gene_name_str)
        else:
            gene_name = gene_name_str

        # Get insertions for this gene, default to empty list if not found
        # Need to find the matching GeneName key since keys are GeneName objects
        insertions = []
        for gene_key, gene_insertions in aa_insertion_set.aa_insertions.items():
            if str(gene_key) == str(gene_name_str):
                insertions = gene_insertions
                break

        # Use per-gene offset from AASequenceSet
        current_offset = aa_sequence_set.offsets[gene_name_str]

        # Create AlignedGene object
        aligned_gene = AlignedGene(
            gene_name=gene_name,
            sequence=sequence,
            offset=current_offset,
            insertions=insertions,
        )

        aligned_genes.append(aligned_gene)

    return aligned_genes
