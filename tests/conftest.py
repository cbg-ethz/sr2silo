"""
This is a configuration file for pytest containing customizations and fixtures.

In VSCode, Code Coverage is recorded in config.xml. Delete this file to reset reporting.
"""

from __future__ import annotations

import contextlib
import tempfile
from pathlib import Path

import pysam
import pytest

from sr2silo.process import bam_to_sam

# Define test data paths
TEST_DATA_DIR = Path(__file__).parent / "data"

INPUT_BAM_PATH = TEST_DATA_DIR / "REF_aln_trim_subsample.bam"
EXPECTED_SAM_PATH = TEST_DATA_DIR / "REF_aln_trim_subsample_expected.sam"

# TODO - centralize the test data.
LARGE_TEST_DATA_DIR = (
    TEST_DATA_DIR
    / "samples_large"
    / "A1_05_2024_10_08/20241024_2411515907"
    / "alignments"
)
# This is a BAM file contains the subset of reads that have insertions
# for effective testing
INPUT_BAM_INSERTIONS_PATH = LARGE_TEST_DATA_DIR / "REF_aln_trim.bam"
EXPECTED_BAM_INSERTIONS_PATH_inserts = (
    LARGE_TEST_DATA_DIR / "REF_aln_trim_subsample_insertions.fasta"
)
EXPECTED_BAM_INSERTIONS_PATH_cleartext = (
    LARGE_TEST_DATA_DIR / "REF_aln_trim_subsample.fasta"
)

# The following file is a BAM file that contains reads that have insertions
# for effective testing
NUC_ALIGNMENT_BAM = TEST_DATA_DIR / "bam" / "combined.bam"


@pytest.fixture
def bam_data() -> Path:
    """Return a sample BAM file path.

    Complementary to sam_data.
    """
    return INPUT_BAM_PATH


@pytest.fixture
def sam_data() -> Path:
    """Return a sample SAM file path.
    This is the expected SAM data for
    the test data in the INPUT_BAM_PATH."""
    return EXPECTED_SAM_PATH


@pytest.fixture
def sam_with_insert_data() -> dict:
    """Return a sample SAM data as a string with insertions."""

    data_expected = dict()
    data_expected["bam_data_fp"] = INPUT_BAM_INSERTIONS_PATH

    # Convert the BAM file to SAM, and read in the SAM file
    with temp_dir() as tmpdir:
        # Convert the BAM file to SAM
        sam_fp = tmpdir / "sam_data.sam"
        bam_to_sam(INPUT_BAM_INSERTIONS_PATH, sam_fp)
        # read in the SAM file
        sam_content = sam_fp.read_text()
        data_expected["sam_data"] = sam_content

    # Read in the insertions file
    insertions_content = EXPECTED_BAM_INSERTIONS_PATH_inserts.read_text()
    data_expected["insertions"] = insertions_content

    # Read in the cleartext file
    cleartext_content = EXPECTED_BAM_INSERTIONS_PATH_cleartext.read_text()
    data_expected["cleartext"] = cleartext_content

    return data_expected


@pytest.fixture
@contextlib.contextmanager
def temp_dir():
    """Return a temporary directory as a Path object.

    This fixture is designed to be used as a context manager.
    """
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
        """Return self."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def fetch(self):
        """Yield a dummy read."""
        yield DummyRead()


class DummyAlignmentFile:
    """A dummy replacement for pysam.AlignmentFile."""

    def __init__(self, *args, **kwargs):
        """Do nothing."""
        pass

    def __enter__(self):
        """Return a DummyBam object."""
        return DummyBam()

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


@pytest.fixture
def dummy_alignment(monkeypatch):
    """Fixture to replace pysam.AlignmentFile with a dummy class."""
    monkeypatch.setattr(pysam, "AlignmentFile", DummyAlignmentFile)
    return Path("dummy.bam")


@pytest.fixture
def primers():
    """Return the primers file path."""
    return Path("./resources/sars-cov-2/primers/primers.yaml")


@pytest.fixture
def timeline():
    """Return the timeline file path."""
    return Path("./tests/data/samples/timeline_A1_05_2024_10_08.tsv")


@pytest.fixture
def sample():
    """Return the sample bam file path."""
    return Path(INPUT_BAM_INSERTIONS_PATH)


@pytest.fixture
def real_sample_files_import_to_loculus(tmp_path, primers, timeline, sample):
    """Get real sample files from the test data directory for
    `sr2silo import-to-loculus`."""
    return {
        "input_file": sample,
        "timeline_file": timeline,
        "primer_file": primers,
        "output_file": tmp_path / "silo_input.ndjson.zst",
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",
        "reference": "sars-cov-2",
    }
