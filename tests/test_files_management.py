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

    def test_csv_positions(self):
        positions_dict = {"pos_ref":[1,2,3], "pos_2":[4,5,6]}
        output_file = '/tmp/'+next(tempfile._get_candidate_names())
        fm.write_positions_to_csv(positions_dict, output_file)
        aligned_positions_dict = fm.get_aligned_positions_dict_from_compotts_output_file(output_file)
        assert(positions_dict==aligned_positions_dict)


if __name__=='__main__':
    unittest.main()



