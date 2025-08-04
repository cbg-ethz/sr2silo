"""Tests for the lapis module."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from sr2silo.loculus.lapis import LapisClient


class TestReferenceJsonToFasta:
    """Test class for the referenceGenomeToFasta function."""

    @pytest.fixture
    def sample_reference_json(self):
        """Sample reference JSON data for testing."""
        return {
            "nucleotideSequences": [{"name": "main", "sequence": "ATCGATCGATCGATCG"}],
            "genes": [
                {
                    "name": "E",
                    "sequence": "MYSFVSEETGTLIVNSVLLFLAFVVFLLVVPDLLV*",
                },
                {
                    "name": "M",
                    "sequence": "MADSNGTITVEELKKLLEQNWITGGILLVQ*",
                },
                {
                    "name": "ORF3a",
                    "sequence": "MDLFMRIFTIGTVTLKQGEIKDAWLIVGVFQSASKIITLKSVPL*",
                },
            ],
        }

    def testreferenceGenomeToFasta_basic(self, sample_reference_json):
        """Test basic functionality of referenceGenomeToFasta."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            nuc_output = tmp_path / "test_nucleotide.fasta"
            aa_output = tmp_path / "test_amino_acid.fasta"

            # Convert dict to JSON string
            json_string = json.dumps(sample_reference_json)

            # Call the function
            LapisClient.referenceGenomeToFasta(
                reference_json_string=json_string,
                nucleotide_out_fp=nuc_output,
                amino_acid_out_fp=aa_output,
            )

            # Verify files were created
            assert nuc_output.exists()
            assert aa_output.exists()

            # Check nucleotide file content
            nuc_content = nuc_output.read_text()
            assert ">main\n" in nuc_content
            assert "ATCGATCGATCGATCG" in nuc_content

            # Check amino acid file content
            aa_content = aa_output.read_text()
            assert ">E\n" in aa_content
            assert ">M\n" in aa_content
            assert ">ORF3a\n" in aa_content

            # Verify stop codons are removed
            assert "*" not in aa_content

    def test_stop_codon_removal(self, sample_reference_json):
        """Test that stop codons (*) are properly removed from amino acid sequences."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            nuc_output = tmp_path / "test_nucleotide.fasta"
            aa_output = tmp_path / "test_amino_acid.fasta"

            json_string = json.dumps(sample_reference_json)

            LapisClient.referenceGenomeToFasta(
                reference_json_string=json_string,
                nucleotide_out_fp=nuc_output,
                amino_acid_out_fp=aa_output,
            )

            aa_content = aa_output.read_text()

            # Check that sequences end without asterisk
            lines = aa_content.strip().split("\n")
            sequence_lines = [
                line for line in lines if line and not line.startswith(">")
            ]

            for line in sequence_lines:
                assert not line.endswith(
                    "*"
                ), f"Sequence should not end with asterisk: {line}"
                assert "*" not in line, f"Sequence should not contain asterisks: {line}"

    def test_invalid_json(self):
        """Test handling of invalid JSON input."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            nuc_output = tmp_path / "test_nucleotide.fasta"
            aa_output = tmp_path / "test_amino_acid.fasta"

            invalid_json = "{ invalid json }"

            with pytest.raises(json.JSONDecodeError):
                LapisClient.referenceGenomeToFasta(
                    reference_json_string=invalid_json,
                    nucleotide_out_fp=nuc_output,
                    amino_acid_out_fp=aa_output,
                )
