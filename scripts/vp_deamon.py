"""Daemon script that processes new samples from the timeline file
   and stores the processed samples in the result directory."""

from __future__ import annotations

import csv
import json
import logging
import sqlite3
import time
from pathlib import Path
import schedule

from vp_transformer import process_directory  # noqa: F401 # isort:skip

logging.basicConfig(level=logging.INFO)

DATABASE_FILE = "processed_files.db"


def initialize_database():
    conn = sqlite3.connect(DATABASE_FILE)
    cursor = conn.cursor()
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS processed_files (
            sample_id TEXT,
            batch_id TEXT,
            PRIMARY KEY (sample_id, batch_id)
        )
    """
    )
    conn.commit()
    conn.close()


def mark_sample_as_processed(sample_id, batch_id):
    conn = sqlite3.connect(DATABASE_FILE)
    cursor = conn.cursor()
    cursor.execute(
        "INSERT OR IGNORE INTO processed_files (sample_id, batch_id) VALUES (?, ?)",
        (sample_id, batch_id),
    )
    conn.commit()
    conn.close()


def is_sample_processed(sample_id, batch_id):
    conn = sqlite3.connect(DATABASE_FILE)
    cursor = conn.cursor()
    cursor.execute(
        "SELECT 1 FROM processed_files WHERE sample_id = ? AND batch_id = ?",
        (sample_id, batch_id),
    )
    result = cursor.fetchone()
    conn.close()
    return result is not None


def read_timeline(timeline_file: Path):
    with timeline_file.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            yield row[0], row[1]


def construct_file_path(sample_id, batch_id):
    return Path(f"samples/{sample_id}/{batch_id}")


def process_new_samples(
    timeline_file: Path, result_dir: Path, nextclade_reference: str
):
    for sample_id, batch_id in read_timeline(timeline_file):
        if not is_sample_processed(sample_id, batch_id):
            logging.info(f"Processing new sample: {sample_id}, batch: {batch_id}")
            file_path = construct_file_path(sample_id, batch_id)
            process_directory(file_path, result_dir, nextclade_reference, timeline_file)
            mark_sample_as_processed(sample_id, batch_id)


def load_config(config_file: Path) -> dict:
    with config_file.open() as f:
        return json.load(f)


def main():
    # Load the configuration
    config = load_config("vp_config.json")

    timeline_file = Path(config["timeline_file"])
    result_dir = Path(config["result_dir"])
    nextclade_reference = config["nextclade_reference"]

    initialize_database()
    schedule.every(10).minutes.do(
        process_new_samples, timeline_file, result_dir, nextclade_reference
    )

    while True:
        schedule.run_pending()
        time.sleep(1)


if __name__ == "__main__":
    main()
