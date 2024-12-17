"""Testing the S3 module."""

from __future__ import annotations

from sr2silo.s3 import compress_bz2


def test_compress_bz2(tmp_path):
    """Test the compress_bz2 function."""
    # Create a test file
    test_file = tmp_path / "test.txt"
    test_file.write_text(
        """
                        This is a test file – it should compress well.
                        This is a test file – it should compress well.
                        This is a test file – it should compress well.
                         """
    )

    # Compress the test file
    compressed_file = tmp_path / "test.txt.bz2"
    compress_bz2(test_file, compressed_file)

    # Check that the compressed file exists
    assert compressed_file.exists()
    # Check that the compressed file is smaller than the original file
    assert compressed_file.stat().st_size < test_file.stat().st_size
