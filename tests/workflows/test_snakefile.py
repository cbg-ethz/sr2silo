import os
import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path


from common import OutputChecker

def test_get_coverage_for_batch():
    """
    Test the get_coverage_for_batch rule.
    using test data from sars-cov-2.

    This version of the test automatically finds the necessary files.
    """
    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        workdir.mkdir(exist_ok=True)

        # Create necessary subdirectories
        (workdir / "config").mkdir(exist_ok=True)
        (workdir / "data").mkdir(exist_ok=True)
        (workdir / "results").mkdir(exist_ok=True)
        (workdir / "scripts").mkdir(exist_ok=True)  # for the script

        # Define paths
        mock_data_path = Path("workflow/.tests/unit/amplicon_cov/data")
        expected_path = Path("workflow/.tests/unit/amplicon_cov/expected")
        config_path = Path("workflow/.tests/unit/amplicon_cov/amplicon_cov.yaml")
        script_path = Path("scripts/amplicon_covs.py")

        # Copy config to the temporary workdir
        wrk_config_path = workdir / "config" / config_path.name
        shutil.copy(config_path, wrk_config_path)

        # Copy mock data to the temporary workdir
        wrk_mock_data_path = Path(workdir, "data")
        shutil.copytree(mock_data_path, wrk_mock_data_path, dirs_exist_ok=True)
        shutil.copy(script_path, workdir / "scripts" / script_path.name)

        # Print the contents of the workdir

        # Print the contents of the current directory
        print_directory_contents(workdir)
        print_directory_contents(wrk_mock_data_path)

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
                "20200729/cov_heatmap.pdf",
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