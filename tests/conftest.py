"""
This is a configuration file for pytest containing customizations and fixtures.

In VSCode, Code Coverage is recorded in config.xml. Delete this file to reset reporting.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List

import pysam
import pytest
from _pytest.nodes import Item

from sr2silo.process import bam_to_sam


# TODO: Add custom markers here
def pytest_collection_modifyitems(items: List[Item]):
    """Add custom markers to tests."""
    for item in items:
        if "spark" in item.nodeid:
            item.add_marker(pytest.mark.spark)
        elif "_int_" in item.nodeid:
            item.add_marker(pytest.mark.integration)


@pytest.fixture
def unit_test_mocks(monkeypatch: None):
    """Include Mocks here to execute all commands offline and fast."""
    pass


# Define test data paths
TEST_DATA_DIR = Path(__file__).parent / "data"

INPUT_BAM_PATH = TEST_DATA_DIR / "REF_aln_trim_subsample.bam"
EXPECTED_SAM_PATH = TEST_DATA_DIR / "REF_aln_trim_subsample_expected.sam"


LARGE_TEST_DATA_DIR = (
    TEST_DATA_DIR
    / "samples_large"
    / "A1_05_2024_10_08/20241024_2411515907"
    / "alignments"
)
INPUT_BAM_INSERTIONS_PATH = LARGE_TEST_DATA_DIR / "REF_aln_trim.bam"
EXPECTED_BAM_INSERTIONS_PATH_inserts = (
    LARGE_TEST_DATA_DIR / "REF_aln_trim_subsample_insertions.fasta"
)
EXPECTED_BAM_INSERTIONS_PATH_cleartext = (
    LARGE_TEST_DATA_DIR / "REF_aln_trim_subsample.fasta"
)


"""Returns a sample BAM data path and its corresponding SAM data as a string."""


@pytest.fixture
def bam_data() -> dict:
    """Return a sample BAM data path and its corresponding SAM data as a string."""
    dict_data = dict()
    # read in these files as test
    dict_data["bam_path"] = INPUT_BAM_PATH
    with open(EXPECTED_SAM_PATH) as f:
        dict_data["sam_data"] = f.read()

    return dict_data


@pytest.fixture
def sam_data():
    """Return a sample SAM data as a string."""
    return bam_to_sam(INPUT_BAM_PATH)


@pytest.fixture
def sam_with_insert_data():
    """Return a sample SAM data as a string with insertions."""

    data_expected = dict()
    data_expected["bam_data_fp"] = INPUT_BAM_INSERTIONS_PATH
    data_expected["sam_data"] = bam_to_sam(INPUT_BAM_INSERTIONS_PATH)
    # Read in the insertions file
    insertions_content = EXPECTED_BAM_INSERTIONS_PATH_inserts.read_text()
    data_expected["insertions"] = insertions_content

    # Read in the cleartext file
    cleartext_content = EXPECTED_BAM_INSERTIONS_PATH_cleartext.read_text()
    data_expected["cleartext"] = cleartext_content

    return data_expected


@pytest.fixture
def temp_dir():
    """Return a temporary directory as a Path object."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield Path(tmpdirname)


class DummyRead:
    """A dummy read object for testing purposes."""

    def __init__(self):
        self.query_name = "read1"
        self.query_sequence = "ACTG"
        self.query_qualities = [30, 31, 32, 33]
        self.reference_start = 100
        self.cigartuples = [(0, 2), (1, 1), (0, 1)]
        self.is_unmapped = False


class DummyBam:
    """A dummy BAM file that yields one dummy read."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def fetch(self):
        """Yield a dummy read."""
        yield DummyRead()


class DummyAlignmentFile:
    """A dummy replacement for pysam.AlignmentFile."""

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return DummyBam()

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


@pytest.fixture
def dummy_alignment(monkeypatch):
    """Fixture to replace pysam.AlignmentFile with a dummy class."""
    monkeypatch.setattr(pysam, "AlignmentFile", DummyAlignmentFile)
    return Path("dummy.bam")
