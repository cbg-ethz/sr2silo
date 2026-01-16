"""Test the process_sample rule with RSV-A data and reference filtering."""

from __future__ import annotations

import os
import shutil
import subprocess as sp
import sys
import time
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))


def test_process_sample_rsva_with_reference_filtering():
    """Test the process_sample rule with RSV-A data and reference filtering.

    This test verifies that:
    1. The snakemake workflow can process RSV-A samples
    2. The REFERENCE_ACCESSION config parameter is passed correctly
    3. Reads are filtered to only include those aligned to the target reference

    The test uses H3_16_2025_11_15 sample which has reads aligned to both:
    - EPI_ISL_412866 (RSV-A): 3086 reads
    - EPI_ISL_1653999 (RSV-B): 48711 reads

    With REFERENCE_ACCESSION: "EPI_ISL_412866", only RSV-A reads should be processed.
    """
    print("Starting RSV-A process_sample test with reference filtering")
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = Path("tests/snakemake/process_sample_rsva/data")
        config_path = Path("tests/snakemake/process_sample_rsva/config.yaml")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Copy the config file to the workdir.
        shutil.copytree(config_path.parent, workdir / "workflow")

        # Copy the main Snakefile if it exists
        main_snakefile = Path("workflow/Snakefile")
        if main_snakefile.exists():
            shutil.copy2(main_snakefile, workdir / "workflow" / "Snakefile")

        # Copy the conda environment files
        envs_dir = Path("workflow/envs")
        if envs_dir.exists():
            shutil.copytree(envs_dir, workdir / "workflow" / "envs")

        # Copy the "resources" folder to the workdir.
        shutil.copytree("resources", workdir / "resources")

        # make dir for results
        os.makedirs(workdir / "results", exist_ok=True)

        # Run the test job with a timeout
        command = [
            "python",
            "-m",
            "snakemake",
            "results/sampleId-H3_16_2025_11_15.ndjson.zst",
            "-f",
            "-j1",
            "--directory",
            str(workdir),
        ]

        print(f"Running command: {' '.join(command)}")
        print(f"Working directory: {workdir}")

        # Set CI environment variable for the subprocess
        env = os.environ.copy()
        env["CI"] = "true"
        try:
            # Add timeout to prevent indefinite hanging (5 minutes)
            process = sp.Popen(
                command, stdout=sp.PIPE, stderr=sp.PIPE, text=True, env=env
            )

            # Wait for process with timeout
            timeout = 300  # 5 minutes
            start_time = time.time()

            while process.poll() is None:
                if time.time() - start_time > timeout:
                    process.kill()
                    raise TimeoutError(f"Process timed out after {timeout} seconds")
                time.sleep(1)

            stdout, stderr = process.communicate()

            # Check if process completed successfully
            if process.returncode != 0:
                print(f"Process failed with return code {process.returncode}")
                print(f"STDOUT: {stdout}")
                print(f"STDERR: {stderr}")

                # Check for log file
                log_file = (
                    workdir / "logs/sr2silo/process_sample/"
                    "sampleId_H3_16_2025_11_15.log"
                )
                if log_file.exists():
                    with open(log_file) as f:
                        log_content = f.read()
                        print(f"Log file content:\n{log_content}")
                else:
                    print(f"Log file {log_file} does not exist.")

                raise AssertionError(
                    f"Snakemake process failed with return code {process.returncode}"
                )

            # Verify the reference filtering was applied by checking log output
            # The log should contain filtering statistics
            log_file = (
                workdir / "logs/sr2silo/process_sample/sampleId_H3_16_2025_11_15.log"
            )
            if log_file.exists():
                with open(log_file) as f:
                    log_content = f.read()
                    print(f"Log file content:\n{log_content}")

                    # Verify reference filtering was applied
                    assert "Reference filtering for 'EPI_ISL_412866'" in log_content, (
                        "Expected reference filtering log message not found"
                    )
                    # Should have filtering statistics with percentages
                    # Note: exact read counts depend on paired-end merging
                    # After merging, RSV-A reads are ~6% and RSV-B are ~94%
                    assert "reads kept" in log_content, (
                        "Expected 'reads kept' in filtering message"
                    )
                    assert "filtered out" in log_content, (
                        "Expected 'filtered out' in filtering message"
                    )
                    # Check that most reads were filtered (RSV-B predominates)
                    assert "94.0%" in log_content or "94." in log_content, (
                        "Expected ~94% filtering ratio (RSV-B reads filtered)"
                    )

        except sp.CalledProcessError as e:
            print(f"CalledProcessError: {e}")
            log_file = (
                workdir / "logs/sr2silo/process_sample/sampleId_H3_16_2025_11_15.log"
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

        # Verify output file was created
        output_file = workdir / "results/sampleId-H3_16_2025_11_15.ndjson.zst"
        if not output_file.exists():
            print(f"Output file {output_file} was not created!")
            # List what files were created in the results directory
            print("Files in results directory:")
            for file in (workdir / "results").glob("*"):
                print(f"  {file.name}")
            raise FileNotFoundError(f"Expected output file not found: {output_file}")

        # Verify output file has content (not empty)
        assert output_file.stat().st_size > 0, "Output file is empty"

        print("RSV-A process_sample test with reference filtering PASSED")
