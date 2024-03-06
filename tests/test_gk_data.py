# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import unittest
import tempfile
import datetime
import calendar
import os
from genome_kit import gk_data
from . import test_data_dir


class TestGkData(unittest.TestCase):
    def test_get_file(self):
        assert not os.path.exists(
            os.path.join(test_data_dir, "gencode.vM30.cfg")
        )

        # Download a test file and check if the content is correct
        path = gk_data.get_file('gencode.vM30.cfg')
        try:
            self.assertTrue(os.path.isfile(path))
            with open(path) as f:
                content = f.read().rstrip()
            self.assertEqual(content, "refg=mm39")

            # Get the test file again, when it already exists
            path = gk_data.get_file('gencode.vM30.cfg')
        finally:
            os.remove(path)


class TestWGet(unittest.TestCase):
    """Test the `gk_data.wget` function."""

    @classmethod
    def setUpClass(cls):
        # Two public-facing small files, with url1 SMALLER and NEWER than url2
        cls.url1 = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10090/108/GCF_000001635.26_GRCm38.p6/annotation_hashes.txt"
        cls.url2 = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10090/108/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_feature_count.txt.gz"
        cls.url1_ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10090/108/GCF_000001635.26_GRCm38.p6/annotation_hashes.txt"
        cls.url2_ftp = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/10090/108/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_feature_count.txt.gz"
        cls.file1 = cls.url1.split("/")[-1]
        cls.file2 = cls.url2.split("/")[-1]
        cls.time1 = calendar.timegm(datetime.datetime(2019, 8, 10, 14, 59, 57).timetuple())  # Newer UTC time
        cls.time2 = calendar.timegm(datetime.datetime(2019, 8, 10, 14, 58, 40).timetuple())  # Older UTC time
        cls.size1 = 411  # Smaller file size
        cls.size2 = 1287  # Larger file size
        cls.tmpdir = tempfile.mktemp()  # Deliberately mktemp (not mkdtemp) to avoid creating dir
        cls.tmpfile = tempfile.mktemp()

    def tearDown(self):
        # Each test starts with tmpfile and non-existent
        if os.path.exists(self.tmpfile): os.unlink(self.tmpfile)
        if os.path.exists(self.tmpdir): os.unlink(self.tmpdir)

    def test_default(self):
        # Check that dst defaults to temp dir, not local dir
        dst = gk_data.wget(self.url1)
        self.assertTrue(dst.startswith(tempfile.gettempdir()))
        self.assertTrue(os.path.isfile(dst))
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)
        os.unlink(dst)

        # Check that specific dst works
        dst = gk_data.wget(self.url1, self.tmpfile)
        self.assertEqual(dst, self.tmpfile)
        self.assertTrue(os.path.isfile(dst))
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)
        os.unlink(dst)

        # Check that creating a directory works
        tmpfile = os.path.join(self.tmpdir, "test.data")
        dst = gk_data.wget(self.url1, tmpfile)
        self.assertEqual(dst, tmpfile)
        self.assertTrue(os.path.isfile(dst))
        self.assertTrue(os.path.isdir(self.tmpdir))
        os.unlink(dst)
        os.rmdir(self.tmpdir)

        # Check that downloading overtop of a file works
        dst = gk_data.wget(self.url1, self.tmpfile)
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)
        dst = gk_data.wget(self.url2, self.tmpfile)
        self.assertEqual(os.path.getsize(dst), self.size2)
        self.assertEqual(os.path.getmtime(dst), self.time2)
        os.unlink(dst)

    def test_timestamping(self):
        # Check that url1 (newer) replaces url2 (older), but not vice versa
        dst = gk_data.wget(self.url2, self.tmpfile)
        self.assertEqual(os.path.getsize(dst), self.size2)
        self.assertEqual(os.path.getmtime(dst), self.time2)
        dst = gk_data.wget(self.url1, self.tmpfile, timestamping=True)
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)
        dst = gk_data.wget(self.url2, self.tmpfile, timestamping=True)
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)

        # Test FTP server where timestamps are not returned as part of the header
        dst = gk_data.wget(self.url2_ftp, self.tmpfile)
        self.assertEqual(os.path.getsize(dst), self.size2)
        self.assertEqual(os.path.getmtime(dst), self.time2)
        dst = gk_data.wget(self.url1_ftp, self.tmpfile, timestamping=True)
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)
        dst = gk_data.wget(self.url2_ftp, self.tmpfile, timestamping=True)
        self.assertEqual(os.path.getsize(dst), self.size1)
        self.assertEqual(os.path.getmtime(dst), self.time1)

    def test_try_ftp_last_modified_date(self):
        # Date should match
        self.assertEqual(gk_data._try_ftp_last_modified_date(self.url1_ftp), self.time1)

        # Getting date of non-FTP URL should return None
        self.assertIs(gk_data._try_ftp_last_modified_date(self.url1), None)

        # Getting modified date of a directory is not implemented
        with self.assertRaises(BaseException):
            gk_data._try_ftp_last_modified_date("ftp://hgdownload.cse.ucsc.edu/goldenPath/")

    def test_no_filesize(self):
        # Check that url with no Content-length attribute still works
        # (Unfortunately can't test this with progressbar, currently.)
        dst = gk_data.wget("http://www.google.com/favicon.ico", self.tmpfile)
        self.assertTrue(os.path.isfile(dst))

    def test_errors(self):
        # Invalid arg types
        with self.assertRaises(TypeError):
            gk_data.wget(242)
        with self.assertRaises(TypeError):
            gk_data.wget(self.url1, 5324)

        # Invalid URL
        with self.assertRaises(IOError):
            gk_data.wget("http://www.invalid-domain-name-123-xyz.com/", self.tmpfile)
        self.assertFalse(os.path.exists(self.tmpfile))
