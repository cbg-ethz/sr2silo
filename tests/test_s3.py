import boto3
import pytest
from moto import mock_aws
from sr2silo.s3 import upload_file_to_s3, download_file_from_s3, list_files_in_bucket

bucket_name = "test-bucket"

@pytest.fixture
def setUp(self):
    self.mock_aws = mock_aws()
    self.mock_aws.start()

    # you can use boto3.client("s3") if you prefer
    s3 = boto3.resource("s3")
    s3_client = s3.boto3.client("s3", region_name="us-east-1")
    s3_client.create_bucket(Bucket=bucket_name)


def test_upload_file_to_s3(s3_client, tmp_path):
    file_content = b'This is a test file.'
    file_name = tmp_path / 'test_file.txt'
    file_name.write_bytes(file_content)
    bucket_name = 'my-test-bucket'

    upload_file_to_s3(str(file_name), bucket_name)

    response = s3_client.get_object(Bucket=bucket_name, Key='test_file.txt')
    assert response['Body'].read() == file_content

def test_download_file_from_s3(s3_client, tmp_path):
    file_content = b'This is a test file.'
    file_name = 'test_file.txt'
    bucket_name = 'my-test-bucket'
    s3_client.put_object(Bucket=bucket_name, Key=file_name, Body=file_content)

    download_path = tmp_path / 'downloaded_test_file.txt'
    download_file_from_s3(bucket_name, file_name, str(download_path))

    assert download_path.read_bytes() == file_content

def test_list_files_in_bucket(s3_client):
    file_content = b'This is a test file.'
    file_name = 'test_file.txt'
    bucket_name = 'my-test-bucket'
    s3_client.put_object(Bucket=bucket_name, Key=file_name, Body=file_content)

    files = list_files_in_bucket(bucket_name)

    assert file_name in files