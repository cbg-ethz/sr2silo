"""Tests for the translation module."""


from __future__ import annotations

from pathlib import Path

from sr2silo.translation import translate


def test_translate():
    """Test the translation function."""

    translate(
        [Path("tests/data/merged_expected.fasta")],
        Path("output/"),
        "nextstrain/sars-cov-2/XBB",
    )
    assert True
