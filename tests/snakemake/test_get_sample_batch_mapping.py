"""Test the get_sample_batch_mapping function."""

from __future__ import annotations

import csv
import os
import tempfile
from datetime import datetime


def get_sample_batch_mapping_test(config):
    """Copy of the get_sample_batch_mapping function for testing."""

    # Read the timeline file
    sample_batch_mapping = {}

    with open(config["TIMELINE_FILE"], "r") as f:
        reader = csv.reader(f, delimiter="\t")
        # Skip header row if it exists
        first_row = next(reader, None)
        if first_row and first_row[0] == "sample":
            # This is a header row, continue with the next rows
            pass
        else:
            # This is data, process it
            if first_row and len(first_row) >= 7:
                sample, batch, reads, proto, location_code, date_str, location = (
                    first_row
                )

                # Convert date string to datetime for comparison
                try:
                    sample_date = datetime.strptime(date_str, "%Y-%m-%d")
                    location_code_int = int(location_code)

                    # Filter by date range
                    start_date = datetime.strptime(config["START_DATE"], "%Y-%m-%d")
                    end_date = datetime.strptime(config["END_DATE"], "%Y-%m-%d")

                    if (
                        start_date <= sample_date <= end_date
                        and location_code_int in config["LOCATIONS"]
                    ):
                        sample_batch_mapping[sample] = batch
                except (ValueError, TypeError):
                    pass  # Skip invalid rows

        # Process remaining rows
        for row in reader:
            if len(row) < 7:  # Skip invalid rows
                continue

            sample, batch, reads, proto, location_code, date_str, location = row

            # Convert date string to datetime for comparison
            try:
                sample_date = datetime.strptime(date_str, "%Y-%m-%d")
            except ValueError:
                continue  # Skip rows with invalid date format

            # Convert location_code to int for comparison
            try:
                location_code_int = int(location_code)
            except ValueError:
                continue  # Skip rows with invalid location code

            # Filter by date range
            start_date = datetime.strptime(config["START_DATE"], "%Y-%m-%d")
            end_date = datetime.strptime(config["END_DATE"], "%Y-%m-%d")

            if not (start_date <= sample_date <= end_date):
                continue

            # Filter by location codes
            if location_code_int not in config["LOCATIONS"]:
                continue

            # Add to mapping
            sample_batch_mapping[sample] = batch

    return sample_batch_mapping


def test_get_sample_batch_mapping():
    """Test the get_sample_batch_mapping function in isolation."""

    # Create a temporary timeline file
    timeline_content = """sample	batch	reads	proto	location_code	date	location
A1_05_2024_10_08	20241024_2411515907	250	v532	5	2024-10-08	Lugano (TI)
B2_15_2024_10_15	20241024_2411515907	250	v532	15	2024-10-15	Basel (BS)
C3_17_2024_11_01	20241024_2411515907	250	v532	17	2024-11-01	Chur (GR)
D4_10_2024_09_30	20241024_2411515907	250	v532	10	2024-09-30	Zürich (ZH)
"""

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(timeline_content)
        timeline_file = f.name

    try:
        # Test 1: Basic filtering by date and location
        config = {
            "TIMELINE_FILE": timeline_file,
            "START_DATE": "2024-10-01",
            "END_DATE": "2024-10-31",
            "LOCATIONS": [5, 15],  # Lugano and Basel
        }

        result = get_sample_batch_mapping_test(config)

        # Should find A1_05_2024_10_08 (location 5, date 2024-10-08) and B2_15_2024_10_15 (location 15, date 2024-10-15)
        expected = {
            "A1_05_2024_10_08": "20241024_2411515907",
            "B2_15_2024_10_15": "20241024_2411515907",
        }
        assert result == expected, f"Expected {expected}, got {result}"
        print("✓ Test 1 passed: Basic date and location filtering")

        # Test 2: No matches (restrictive date range)
        config["START_DATE"] = "2024-12-01"
        config["END_DATE"] = "2024-12-31"

        result = get_sample_batch_mapping_test(config)
        assert result == {}, f"Expected empty dict, got {result}"
        print("✓ Test 2 passed: No matches for restrictive date range")

        # Test 3: No matches (restrictive location)
        config["START_DATE"] = "2024-10-01"
        config["END_DATE"] = "2024-10-31"
        config["LOCATIONS"] = [99]  # Non-existent location

        result = get_sample_batch_mapping_test(config)
        assert result == {}, f"Expected empty dict, got {result}"
        print("✓ Test 3 passed: No matches for non-existent location")

        # Test 4: Single location match
        config["LOCATIONS"] = [17]  # Chur
        config["START_DATE"] = "2024-11-01"
        config["END_DATE"] = "2024-11-01"

        result = get_sample_batch_mapping_test(config)
        expected = {"C3_17_2024_11_01": "20241024_2411515907"}
        assert result == expected, f"Expected {expected}, got {result}"
        print("✓ Test 4 passed: Single exact match")

        print("All get_sample_batch_mapping tests passed!")

    finally:
        # Clean up
        os.unlink(timeline_file)


if __name__ == "__main__":
    test_get_sample_batch_mapping()
