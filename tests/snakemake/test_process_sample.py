"""Test the process_sample rule."""

from __future__ import annotations

import os
import shutil
import subprocess as sp
import sys
import time
from pathlib import Path
from tempfile import TemporaryDirectory

import zstandard as zstd

sys.path.insert(0, os.path.dirname(__file__))


# Nota bene: The output tested against is not validated, yet a failure in test
# notes a change of output here.
def test_process_sample():
    """Test the process_sample rule."""

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"  # Ensure Path object
        data_path = Path("tests/snakemake/process_sample/data")
        expected_path = Path(
            "tests/snakemake/process_sample/expected/results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst"
        )
        config_path = Path("tests/snakemake/process_sample/config.yaml")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy the config file to the workdir.
        shutil.copytree(config_path.parent, workdir / "workflow")

        # Copy the "resources" folder to the workdir.
        shutil.copytree("resources", workdir / "resources")

        # make dir for results
        os.makedirs(workdir / "results", exist_ok=True)

        # Run the test job with a timeout
        command = [
            "python",
            "-m",
            "snakemake",
            "results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst",
            "-f",
            "-j1",
            "--target-files-omit-workdir-adjustment",
            "--directory",
            str(workdir),  # Convert Path to string
        ]

        print(f"Running command: {' '.join(command)}")
        print(f"Working directory: {workdir}")

        try:
            # Add timeout to prevent indefinite hanging (5 minutes)
            process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, text=True)

            # Wait for process with timeout
            timeout = 300  # 5 minutes
            start_time = time.time()

            while process.poll() is None:
                if time.time() - start_time > timeout:
                    process.kill()
                    raise TimeoutError(f"Process timed out after {timeout} seconds")
                time.sleep(1)

            # Check if process completed successfully
            if process.returncode != 0:
                stdout, stderr = process.communicate()
                print(f"Process failed with return code {process.returncode}")
                print(f"STDOUT: {stdout}")
                print(f"STDERR: {stderr}")

                # Check for log file
                log_file = (
                    workdir / "logs/sr2silo/process_sample/"
                    "sampleId_A1_05_2024_10_08_batchId_20241024_2411515907.log"
                )
                if log_file.exists():
                    with open(log_file) as f:
                        log_content = f.read()
                        print(f"Log file content:\n{log_content}")
                else:
                    print(f"Log file {log_file} does not exist.")

                raise sp.CalledProcessError(process.returncode, command)

        except sp.CalledProcessError as e:
            log_file = (
                workdir / "logs/sr2silo/process_sample/"
                "sampleId_A1_05_2024_10_08_batchId_20241024_2411515907.log"
            )
            if log_file.exists():
                with open(log_file) as f:
                    print(f.read())
            else:
                print(f"Log file {log_file} does not exist.", file=sys.stderr)

            # Check if the Snakefile exists
            snakefile = workdir / "workflow" / "Snakefile"
            if not snakefile.exists():
                print(f"Snakefile not found at {snakefile}", file=sys.stderr)

            # List directory contents to help diagnose issues
            print(f"Contents of {workdir}:")
            for path in workdir.glob("**/*"):
                if path.is_file():
                    print(f"  {path.relative_to(workdir)}")

            raise e

        # The rest of the function remains the same
        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        # common.OutputChecker(data_path, expected_path, workdir).check()

        output_file = (
            workdir
            / "results/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.ndjson.zst"
        )
        if not output_file.exists():
            print(f"Output file {output_file} was not created!")
            # List what files were created in the results directory
            print("Files in results directory:")
            for file in (workdir / "results").glob("*"):
                print(f"  {file.name}")
            raise FileNotFoundError(f"Expected output file not found: {output_file}")

        # Decompress both files and compare them as strings as they are not binary
        with open(output_file, "rb") as f:
            decompressor = zstd.ZstdDecompressor()
            with decompressor.stream_reader(f) as reader:
                generated_content = reader.read().decode("utf-8")

        with open(expected_path, "rb") as f:
            decompressor_expected = zstd.ZstdDecompressor()
            with decompressor_expected.stream_reader(f) as reader:
                expected_content = reader.read().decode("utf-8")

        error_message = (
            "The contents of the generated and expected files do not match.\n\n"
            "Generated content:\n"
            f"{generated_content}\n\n"
            "Expected content:\n"
            f"{expected_content}"
        )

        assert generated_content == expected_content, error_message
