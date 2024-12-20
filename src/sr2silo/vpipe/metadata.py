"""Extract metadata from V-Pipe Filenaming Conventions."""

from __future__ import annotations


def sample_id_decoder(sample_id: str) -> dict:
    """Decode the sample ID into individual components.

    Args:
        sample_id (str): The sample ID to decode.

    Returns:
        dict: A dictionary containing the decoded components.
              containing the following keys:
                - sequencing_well_position (str : sequencing well position)
                - location_code (int : code of the location)
                - sampling_date (str : date of the sampling)
    """
    components = sample_id.split("_")
    # Assign components to meaningful variable names
    well_position = components[0]  # A1
    location_code = components[1]  # 10
    sampling_date = f"{components[2]}-{components[3]}-{components[4]}"  # 2024-09-30
    return {
        "sequencing_well_position": well_position,
        "location_code": location_code,
        "sampling_date": sampling_date,
    }
