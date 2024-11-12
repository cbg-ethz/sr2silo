"""Tests for the translation module."""


from __future__ import annotations

from sr2silo.translation import translate


def test_translate():
    """Test the translation function."""
    translate("tests/data/merged_expected.fasta", "output/", "data/sars-cov-2")
    assert True
