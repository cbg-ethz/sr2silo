"""Test the process_sample rule."""

from __future__ import annotations

import os
import shutil
import subprocess as sp
import sys
import time
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))


# Nota bene: The output tested against is not validated, yet a failure in test
# notes a change of output here.
def test_process_sample():
    """Test the process_sample rule. - does not validate the output."""
    print("starting")
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"  # Ensure Path object
        data_path = Path("tests/snakemake/process_sample/data")
        config_path = Path("tests/snakemake/process_sample/config.yaml")

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
            "results/sampleId-A1_05_2024_10_08.ndjson.zst",
            "-f",
            "-j1",
            "--use-conda",  # Force use of conda environments
            "--conda-frontend",
            "conda",  # Use conda instead of mamba
            "--target-files-omit-workdir-adjustment",
            "--directory",
            str(workdir),  # Convert Path to string
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

            # Check if process completed successfully
            if process.returncode != 0:
                stdout, stderr = process.communicate()
                print(f"Process failed with return code {process.returncode}")
                print(f"STDOUT: {stdout}")
                print(f"STDERR: {stderr}")

                # Check for log file
                log_file = (
                    workdir / "logs/sr2silo/process_sample/"
                    "sampleId_A1_05_2024_10_08.log"
                )
                if log_file.exists():
                    with open(log_file) as f:
                        log_content = f.read()
                        print(f"Log file content:\n{log_content}")
                else:
                    print(f"Log file {log_file} does not exist.")

        except sp.CalledProcessError as e:
            print(f"CalledProcessError: {e}")
            log_file = (
                workdir / "logs/sr2silo/process_sample/" "sampleId_A1_05_2024_10_08.log"
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

        output_file = workdir / "results/sampleId-A1_05_2024_10_08.ndjson.zst"
        if not output_file.exists():
            print(f"Output file {output_file} was not created!")
            # List what files were created in the results directory
            print("Files in results directory:")
            for file in (workdir / "results").glob("*"):
                print(f"  {file.name}")
            raise FileNotFoundError(f"Expected output file not found: {output_file}")
