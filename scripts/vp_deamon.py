"""Daemon script that processes new samples from the timeline file
   and stores the processed samples in the result directory.
   """

from __future__ import annotations

import csv
import datetime
import json
import logging
import shutil
import sqlite3
import time
from pathlib import Path

import click
import schedule

from vp_transformer import process_directory  # noqa: F401 # isort:skip

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)

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


def mark_sample_as_processed(database_file: Path, sample_id, batch_id):
    conn = sqlite3.connect(database_file)
    cursor = conn.cursor()
    cursor.execute(
        "INSERT OR IGNORE INTO processed_files (sample_id, batch_id) VALUES (?, ?)",
        (sample_id, batch_id),
    )
    conn.commit()
    conn.close()


def is_sample_processed(database_file: Path, sample_id, batch_id):
    conn = sqlite3.connect(database_file)
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


def construct_file_path(sample_dir: Path, sample_id: str, batch_id: str) -> Path:
    return sample_dir / sample_id / batch_id


def process_new_samples(
    database_file: Path,
    sample_dir: Path,
    timeline_file: Path,
    result_dir: Path,
    nextclade_reference: str,
):
    for sample_id, batch_id in read_timeline(timeline_file):
        if not is_sample_processed(database_file, sample_id, batch_id):
            logging.info(f"Processing new sample: {sample_id}, batch: {batch_id}")
            file_path = construct_file_path(sample_dir, sample_id, batch_id)
            output_dir = result_dir / sample_id / batch_id
            output_dir.mkdir(parents=True, exist_ok=True)
            try:
                process_directory(file_path, output_dir, nextclade_reference, timeline_file)
                mark_sample_as_processed(database_file, sample_id, batch_id)
            except Exception as e:
                logging.error(f"Error processing sample {sample_id}, batch {batch_id}: {e}")


def load_config(config_file: Path) -> dict:
    with config_file.open() as f:
        return json.load(f)


def backup_database(database_file: Path, backup_dir: Path):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file = backup_dir / f"processed_files_{timestamp}.db"
    backup_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(database_file, backup_file)
    logging.info(f"Database backed up to: {backup_file}")


def main():
    # Load the configuration
    logging.info("Loading configuration...")
    config = load_config(Path("scripts/vp_config.json"))

    sample_dir = Path(config["sample_dir"])
    timeline_file = Path(config["timeline_file"])
    result_dir = Path(config["result_dir"])
    nextclade_reference = config["nextclade_reference"]
    database_file = Path(config["database_file"])
    backup_dir = Path(config["backup_dir"])

    logging.info("Initializing database...")
    initialize_database()

    logging.info("Scheduling the sample processing job...")
    schedule.every(1).minutes.do(
        process_new_samples,
        database_file,
        sample_dir,
        timeline_file,
        result_dir,
        nextclade_reference,
    )

    logging.info("Scheduling the database backup job...")
    schedule.every().day.at("02:00").do(backup_database, database_file, backup_dir)

    logging.info("Starting the scheduler...")
    while True:
        schedule.run_pending()
        logging.debug("Waiting for the next scheduled task...")
        time.sleep(10)


if __name__ == "__main__":
    main()
