from __future__ import annotations

import os
import shutil
import subprocess as sp
from pathlib import Path
from tempfile import TemporaryDirectory
import sys

# Add the project root directory to the Python path
sys.path.append(str(Path(__file__).resolve().parents[3]))

from tests.workflow.common import OutputChecker


def test_process_sample():
    """
    Test the process_sample rule.

    This version of the test automatically finds the necessary files.
    """
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        workdir.mkdir(exist_ok=True)

        # Create necessary subdirectories
        (workdir / "config").mkdir(exist_ok=True)
        (workdir / "data").mkdir(exist_ok=True)
        (workdir / "results").mkdir(exist_ok=True)

        # Define paths
        mock_data_path = Path("tests/data/samples_large/A1_05_2024_10_08/20241024_2411515907/alignments/REF_aln_trim.bam")
        expected_path = Path("tests/workflow/process_sample/expected/results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst")
        config_path = Path("tests/workflow/process_sample/config.yaml")

        # Copy config to the temporary workdir
        wrk_config_path = workdir / "config" / config_path.name
        shutil.copy(config_path, wrk_config_path)

        # Copy mock data to the temporary workdir
        wrk_mock_data_path = Path(workdir, "data")
        shutil.copytree(mock_data_path, wrk_mock_data_path, dirs_exist_ok=True)

        # Run the test job

        sp.check_output(
            [
                "snakemake",
                "--snakefile",
                "workflow/rules/amplicon_cov.smk",
                "--configfile",
                str(wrk_config_path),
                "--config",
                "--directory",
                str(workdir),
                "--cores",
                "1",
                "tests/workflow/process_sample/expected/results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst",
            ]
        )

        # Check the output
        # assert (workdir / "results/").exists()

        # show me the full tree of files in the workdir
        for root, dirs, files in os.walk(workdir):
            print(root)
            for file in files:
                print(f"  {file}")

        # Compare output with expected result using the OutputChecker
        checker = OutputChecker(
            workdir / "data",
            expected_path,
            workdir,
            configdir=workdir / "config",
            tolerance=1e-4,
            ignore_files=["pdf", "py", "log"],
        )

        checker.check()
