"""Test the process_sample rule with dry-run first."""

from __future__ import annotations

import os
import shutil
import subprocess as sp
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))


def test_process_sample_dry_run():
    """Test the process_sample rule with dry-run to check syntax."""

    print("Starting process_sample dry-run test")
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = Path("tests/snakemake/process_sample/data")
        config_path = Path("tests/snakemake/process_sample/config.yaml")

        # Copy data to the temporary workdir
        shutil.copytree(data_path, workdir)

        # Copy the config file to the workdir
        shutil.copytree(config_path.parent, workdir / "workflow")

        # Copy the main Snakefile
        main_snakefile = Path("workflow/Snakefile")
        if main_snakefile.exists():
            shutil.copy2(main_snakefile, workdir / "workflow" / "Snakefile")

        # Copy the resources folder to the workdir
        if Path("resources").exists():
            shutil.copytree("resources", workdir / "resources")

        # Make dir for results
        os.makedirs(workdir / "results", exist_ok=True)

        # Run dry-run test
        command = [
            "python",
            "-m",
            "snakemake",
            "results/sampleId-A1_05_2024_10_08.ndjson.zst",
            "-n",  # dry-run
            "-j1",
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
                timeout=60,  # 1 minute timeout for dry-run
                check=True,
            )
            print("Dry-run completed successfully!")
            print(f"STDOUT: {result.stdout}")

            # Check that the process_sample rule was recognized
            assert "process_sample" in result.stdout
            print("✓ process_sample rule syntax is valid")

        except sp.CalledProcessError as e:
            print(f"Command failed with return code {e.returncode}")
            print("STDOUT:", e.stdout)
            print("STDERR:", e.stderr)
            raise e

        except sp.TimeoutExpired:
            print("Command timed out after 60 seconds")
            raise

        print("process_sample dry-run test completed successfully")


if __name__ == "__main__":
    test_process_sample_dry_run()
