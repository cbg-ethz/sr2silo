"""Test the submit_to_loculus rule."""

from __future__ import annotations

import os
import shutil
import subprocess as sp
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))


def test_submit_to_loculus():
    """Test the submit_to_loculus rule."""

    print("Starting submit_to_loculus test")
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"  # Ensure Path object
        data_path = Path("tests/snakemake/submit_to_loculus/data")
        config_path = Path("tests/snakemake/submit_to_loculus/config.yaml")

        # Copy data to the temporary workdir
        shutil.copytree(data_path, workdir)

        # Copy the config file to the workdir
        shutil.copytree(config_path.parent, workdir / "workflow")

        # Copy the Snakefile from the main workflow directory
        main_snakefile = Path("workflow/Snakefile")
        if main_snakefile.exists():
            shutil.copy2(main_snakefile, workdir / "workflow" / "Snakefile")

        # Copy the conda environment files
        envs_dir = Path("workflow/envs")
        if envs_dir.exists():
            shutil.copytree(envs_dir, workdir / "workflow" / "envs")

        # Copy the resources folder to the workdir (if needed)
        if Path("resources").exists():
            shutil.copytree("resources", workdir / "resources")

        # Make dir for uploads output
        os.makedirs(workdir / "results" / "uploads", exist_ok=True)

        # Run the test job with a timeout - using dry-run to avoid S3 upload
        command = [
            "python",
            "-m",
            "snakemake",
            "results/uploads/sampleId-A1_05_2024_10_08_batchId-20241024_2411515907.uploaded",
            "-n",  # dry-run mode - this will check the rule without executing
            "-j1",
            "--use-conda",  # Force use of conda environments
            "--conda-frontend",
            "conda",  # Use conda instead of mamba
            "--target-files-omit-workdir-adjustment",
            "--directory",
            str(workdir),
        ]

        print(f"Running command: {' '.join(command)}")
        print(f"Working directory: {workdir}")

        try:
            result = sp.run(
                command,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minutes timeout
                check=True,
            )
            print("Command executed successfully!")
            print(f"STDOUT: {result.stdout}")

            # Check that the rule was properly parsed and would be executed
            # In dry-run mode,
            # the output file won't be created but the rule should be validated
            assert (
                "submit_to_loculus" in result.stdout
                or "submit_to_loculus" in result.stderr
            )

        except sp.CalledProcessError as e:
            print(f"Command failed with return code {e.returncode}")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr)

            # Check for log file
            log_file = (
                workdir / "logs/sr2silo/submit_to_loculus/"
                "sampleId_A1_05_2024_10_08_batchId_20241024_2411515907.log"
            )
            if log_file.exists():
                with open(log_file) as f:
                    log_content = f.read()
                    print(f"Log file content:\n{log_content}")

            # For a dry-run, we expect it to succeed in parsing the rule
            # but might fail if there are syntax errors
            raise e

        except sp.TimeoutExpired:
            print("Command timed out after 5 minutes")
            raise

        print("submit_to_loculus test completed successfully - rule syntax is valid")


if __name__ == "__main__":
    test_submit_to_loculus()
