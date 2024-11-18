# sr2silo Scripts

This directory contains scripts for processing and managing sample data. Below is an explanation of the two main scripts in this directory.

## vp_deamon.py

`vp_deamon.py` is a daemon script that processes new samples from the timeline file and stores the processed samples in the result directory. It performs the following tasks:

1. **Load Configuration**: Loads configuration settings from a JSON file using Pydantic for validation.
2. **Initialize Database**: Initializes a SQLite database to keep track of processed samples.
3. **Process New Samples**: Reads the timeline file to identify new samples and processes them using the `vp_transformer` module.
4. **Backup Database**: Backs up the database daily to a specified backup directory.
5. **Schedule Tasks**: Uses the `schedule` library to run the sample processing and database backup tasks at specified intervals.

### Usage

To run the daemon script, execute the following command:

```sh
python vp_deamon.py --config scripts/vp_config.json
```
Ensure that the configuration file vp_config.json is present in the scripts directory with the necessary settings.

## vp_transformer.py
`vp_transformer.py` is a script that contains the core processing logic for transforming sample data. This script is used by `vp_deamon.py` to process new samples.

## Usage
This script is typically not run directly. Instead, it is imported and used by `vp_deamon.py`.

## Legacy Notice
The core processing logic in these scripts is based on the dgivec scripts, which were the foundation of this package. These scripts are retained here for legacy reasons and to ensure compatibility with existing workflows.

## Configuration
The configuration file `vp_config.json` should have the following structure:

```json
{
    "sample_dir": "The directory where the samples are stored.",
    "timeline_file": "The path to the timeline file.",
    "result_dir": "The directory where the results are stored.",
    "nextclade_reference": "The reference to use for nextclade.",
    "database_file": "The path to the database file.",
    "backup_dir": "The directory where the backups are stored.",
    "deamon_interval_m": "The interval in minutes to run the daemon."
}
```
- `sample_dir`: The directory where the samples are stored.
- `timeline_file`: The path to the timeline file.
- `result_dir`: The directory where the results are stored.
- `nextclade_reference`: The reference to use for nextclade.
- `database_file`: The path to the database file.
- `backup_dir`: The directory where the backups are stored.
- `deamon_interval_m`: The interval in minutes to run the daemon.

Ensure that all paths are correctly set in the configuration file before running the scripts.
