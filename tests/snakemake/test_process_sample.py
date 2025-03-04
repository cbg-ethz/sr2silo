from __future__ import annotations

import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

import zstandard as zstd

sys.path.insert(0, os.path.dirname(__file__))


def test_process_sample():

    with TemporaryDirectory() as tmpdir:

        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/snakemake/process_sample/data")
        expected_path = PurePosixPath(
            "tests/snakemake/process_sample/expected/results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst"
        )
        config_path = PurePosixPath("tests/snakemake/process_sample/config.yaml")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy the config file to the workdir.
        shutil.copytree(config_path.parent, workdir / "workflow")

        # Copy the "resources" folder to the workdir.
        shutil.copytree("resources", workdir / "resources")

        # make dir for results
        os.makedirs(workdir / "results")

        # dbg
        print(
            "results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst",
                "-f",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--directory",
                workdir,
            ]
        )

        # print logs
        print("logs:")
        # logs/sr2silo/process_sample/sampleId_A1_05_2024_10_08_batchId_20241024_2411515907.log
        with open(
            workdir / "logs/sr2silo/process_sample/"
            "sampleId_A1_05_2024_10_08_batchId_20241024_2411515907.log"
        ) as f:
            print(f.read())

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        # common.OutputChecker(data_path, expected_path, workdir).check()

        # Decompress both files and compare them as strings as they are not binary
        with open(
            workdir / "results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907."
            "ndjson.zst",
            "rb",
        ) as f:
            decompressor = zstd.ZstdDecompressor()
            generated_content = decompressor.decompress(f.read()).decode("utf-8")

        with open(expected_path, "rb") as f:
            expected_content = decompressor.decompress(f.read()).decode("utf-8")

        assert (
            generated_content == expected_content
        ), "The contents of the generated and expected files do not match."
