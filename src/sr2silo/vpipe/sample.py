"""Implements a sample object for V-Pipe."""

from __future__ import annotations

from pathlib import Path

import sr2silo.vpipe.metadata as metadata


class Sample:
    """A sample object for V-Pipe.

    Args:
        sample_id (str): The sample ID.
        organism (str): The organism identifier (e.g., 'covid', 'rsva').
                       Optional, defaults to 'covid' for backward compatibility.
        config_path (Path | None): Optional path to an external timeline columns
                           YAML file. If provided, used instead of the bundled one.
    """

    def __init__(self, sample_id: str, organism: str = "covid", config_path: Path | None = None) -> None:
        self.sample_id = sample_id
        self.organism = organism
        self.config_path = config_path
        self.metadata: dict[str, str] | None = None
        self.timeline: Path | None = None

    def __str__(self) -> str:
        return f"Sample(sample_id={self.sample_id}, organism={self.organism})"

    def enrich_metadata(self, timeline: Path) -> None:
        """Enrich the sample metadata with additional information.

        Args:
            timeline (Path): The path to the timeline file.
        """
        self.timeline = timeline
        self.set_metadata()

    def set_metadata(self) -> None:
        """Get the metadata for the sample."""
        if not self.timeline:
            raise ValueError("Timeline must be set before calling get_metadata")
        self.metadata = metadata.get_metadata(
            sample_id=self.sample_id,
            timeline=self.timeline,
            organism=self.organism,
            config_path=self.config_path,
        )

    def get_metadata(self) -> dict[str, str]:
        """Get the metadata for the sample."""
        if self.metadata is None:
            raise ValueError("Metadata is not set.")
        return self.metadata
