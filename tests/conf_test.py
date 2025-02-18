from __future__ import annotations

from pathlib import Path

import pysam
import pytest


# Dummy read object
class DummyRead:
    def __init__(self):
        self.query_name = "read1"
        self.query_sequence = "ACTG"
        self.query_qualities = [30, 31, 32, 33]
        self.reference_start = 100
        # CIGAR: match 2, insertion 1, match 1
        self.cigartuples = [(0, 2), (1, 1), (0, 1)]
        self.is_unmapped = False


# Dummy BAM file that yields one dummy read
class DummyBam:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def fetch(self):
        yield DummyRead()


# Dummy replacement for pysam.AlignmentFile
class DummyAlignmentFile:
    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return DummyBam()

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


@pytest.fixture
def dummy_alignment(monkeypatch):
    # Replace pysam.AlignmentFile with our dummy class.
    monkeypatch.setattr(pysam, "AlignmentFile", DummyAlignmentFile)
    # Return a dummy BAM file path (not used by the dummy AlignmentFile).
    return Path("dummy.bam")
