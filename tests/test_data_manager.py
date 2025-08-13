import os
import tempfile
import unittest

from genome_kit import DefaultDataManager, GCSDataManager
from genome_kit.data_manager import GKDataFileNotFoundError
from . import test_data_dir


class TestDefaultDataManager(unittest.TestCase):
    @unittest.skipIf('CI' in os.environ, "can't provide AWS credentials in CI")
    def test_get_file(self):
        assert not os.path.exists(os.path.join(test_data_dir, "gencode.vM30.cfg"))

        dm = DefaultDataManager(test_data_dir)
        # Download a test file and check if the content is correct
        path = dm.get_file("gencode.vM30.cfg")
        try:
            self.assertTrue(os.path.isfile(path))
            with open(path) as f:
                content = f.read().rstrip()
            self.assertEqual(content, "refg=mm39")

            # Get the test file again, when it already exists
            path = dm.get_file('gencode.vM30.cfg')
        finally:
            os.remove(path)

    @unittest.skipIf('CI' in os.environ, "can't provide AWS credentials in CI")
    def test_upload_file(self):
        filename = "unittest_upload_file"
        file_contents = "foo"
        assert not os.path.exists(os.path.join(test_data_dir, filename))

        dm = DefaultDataManager(test_data_dir, require_auth=True)
        with tempfile.NamedTemporaryFile("wt", delete=False) as f:
            filepath = f.name
            f.write(file_contents)
        dm.upload_file(filepath, filename)
        path = dm.get_file(filename)
        try:
            self.assertTrue(os.path.isfile(path))
            with open(path) as f:
                content = f.read()
            self.assertEqual(content, file_contents)
        finally:
            os.remove(path)

        # expect an error if the remote file exists
        with tempfile.NamedTemporaryFile("wt", delete=False) as f:
            filepath = f.name
            f.write(file_contents+"bar") # same name, different content
        with self.assertRaises(ValueError) as context:
            dm.upload_file(filepath, filename)
        self.assertIn(f"File '{filename}' already exists", str(context.exception))

        dm.client.delete_object(Bucket=dm._bucket_name, Key=filename)

    @unittest.skipIf('CI' in os.environ, "can't provide AWS credentials in CI")
    def test_nonexistent_genome(self):
        dm = DefaultDataManager(test_data_dir)
        with self.assertRaises(GKDataFileNotFoundError):
            dm.get_file("nonexistent_genome.cfg")


class TestGCSDataManager(unittest.TestCase):
    @unittest.skipIf('CI' in os.environ, "can't provide GCP credentials in CI")
    def test_get_file(self):
        assert not os.path.exists(os.path.join(test_data_dir, "gencode.vM30.cfg"))

        dm = GCSDataManager(test_data_dir, "genomekit-public-dg")
        # Download a test file and check if the content is correct
        path = dm.get_file("gencode.vM30.cfg")
        try:
            self.assertTrue(os.path.isfile(path))
            with open(path) as f:
                content = f.read().rstrip()
            self.assertEqual(content, "refg=mm39")

            # Get the test file again, when it already exists
            path = dm.get_file('gencode.vM30.cfg')
        finally:
            os.remove(path)

    @unittest.skipIf('CI' in os.environ, "can't provide GCP credentials in CI")
    def test_upload_file(self):
        filename = "unittest_upload_file"
        file_contents = "foo"
        assert not os.path.exists(os.path.join(test_data_dir, filename))

        dm = GCSDataManager(test_data_dir, "genomekit-public-dg")
        with tempfile.NamedTemporaryFile("wt", delete=False) as f:
            filepath = f.name
            f.write(file_contents)
        dm.upload_file(filepath, filename)
        path = dm.get_file(filename)
        try:
            self.assertTrue(os.path.isfile(path))
            with open(path) as f:
                content = f.read()
            self.assertEqual(content, file_contents)
        finally:
            os.remove(path)

        # expect an error if the remote file exists
        with tempfile.NamedTemporaryFile("wt", delete=False) as f:
            filepath = f.name
            f.write(file_contents+"bar") # same name, different content
        with self.assertRaises(ValueError) as context:
            dm.upload_file(filepath, filename)
        self.assertIn(f"File '{filename}' already exists", str(context.exception))

        dm.bucket.delete_blob(filename)

    @unittest.skipIf('CI' in os.environ, "can't provide GCP credentials in CI")
    def test_nonexistent_genome(self):
        dm = GCSDataManager(test_data_dir, "genomekit-public-dg")
        with self.assertRaises(GKDataFileNotFoundError):
            dm.get_file("nonexistent_genome.cfg")
