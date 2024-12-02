"""Testing to run Rust based Silo input transformer."""

from __future__ import annotations

import logging

import pytest

logging.basicConfig(level=logging.INFO)

config_path = "silo_input_transformer/config.yaml"


@pytest.mark.skip(
    reason="silo_input_transformer is not public and added as a submodule currently"
)
def test_transform():
    """Test to run Rust based Silo input transformer, fail for error."""

    import silo_input_transformer

    result = silo_input_transformer.run_with_config(config_path)  # type: ignore
    logging.info(f"Result from Rust: {result}")
    return result


if __name__ == "__main__":
    test_transform()
