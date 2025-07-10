"""Implements a sample object for V-Pipe."""

from __future__ import annotations

from pathlib import Path

import sr2silo.vpipe.metadata as metadata


class Sample:
    """A sample object for V-Pipe.

    Args:
        sample_id (str): The sample ID.
        batch_id (str): The batch ID.
    """

    def __init__(self, sample_id: str, batch_id: str) -> None:
        self.sample_id = sample_id
        self.batch_id = batch_id
        self.metadata: dict[str, str] | None = None
        self.timeline: Path | None = None
        self.primers: Path | None = None

    def __str__(self) -> str:
        return f"Sample(sample_id={self.sample_id}, batch_id={self.batch_id})"

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
            raise ValueError(
                "Timeline must be set before calling get_metadata"
            )
        self.metadata = metadata.get_metadata(
            sample_id=self.sample_id,
            batch_id=self.batch_id,
            timeline=self.timeline,
        )

    def get_metadata(self) -> dict[str, str]:
        """Get the metadata for the sample."""
        if self.metadata is None:
            raise ValueError("Metadata is not set.")
        return self.metadata
