#   ---------------------------------------------------------------------------------
#   Copyright (c) Microsoft Corporation. All rights reserved.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   ---------------------------------------------------------------------------------
"""
This is a configuration file for pytest containing customizations and fixtures.

In VSCode, Code Coverage is recorded in config.xml. Delete this file to reset reporting.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import List

import pytest
from _pytest.nodes import Item

from sr2silo.convert import bam_to_sam


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


@pytest.fixture
def sam_data():
    """Return a sample SAM data as a string."""
    return bam_to_sam(INPUT_BAM_PATH)


@pytest.fixture
def temp_dir():
    """Return a temporary directory as a Path object."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        yield Path(tmpdirname)
