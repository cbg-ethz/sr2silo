"""Testing to run Rust based Silo input transformer."""

from __future__ import annotations

import logging

import silo_input_transformer

logging.basicConfig(level=logging.INFO)

config_path = "silo_input_transformer/config.yaml"


def test_transform():
    """Test to run Rust based Silo input transformer, fail for error."""
    result = silo_input_transformer.run_with_config(config_path)  # type: ignore
    logging.info(f"Result from Rust: {result}")
    return result


if __name__ == "__main__":
    test_transform()
