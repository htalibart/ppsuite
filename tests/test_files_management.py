import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils import files_management as fm


class Test_Files_Management(unittest.TestCase):

    def test_file_exists(self):
        existent_filename = SEQ_1CC8
        fm.check_if_file_ok(existent_filename)

    def test_no_file(self):
        nonexistent_filename = "patate.fasta"
        with self.assertRaises(Exception) as context:
            fm.check_if_file_ok(nonexistent_filename)

if __name__=='__main__':
    unittest.main()



