import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from makepotts.potts_object import *
from ppalign.call_ppalign import *
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
        aligned_positions_dict = fm.get_aligned_positions_dict_from_ppalign_output_file(output_file)
        assert(positions_dict==aligned_positions_dict)


    def test_get_aligned_positions_with_gaps_dict_from_ppalign_output_file(self):
        positions_with_gaps_dict={"pos_ref":[0,'-',1], "pos_2":[0,1,'-']}
        output_file = '/tmp/'+next(tempfile._get_candidate_names())
        fm.write_positions_to_csv(positions_with_gaps_dict, output_file)
        aligned_positions_with_gaps_dict = fm.get_aligned_positions_with_gaps_dict_from_ppalign_output_file(output_file)
        assert(positions_with_gaps_dict==aligned_positions_with_gaps_dict)


    def test_get_seqs_aligned_in_fasta_file(self):
        potts_folder_1 = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        potts_folder_2 = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        obj1 = Potts_Object.from_sequence_alone(potts_folder_1, DUMMY_SEQUENCE_FILE_1, "one_submat")
        obj2 = Potts_Object.from_sequence_alone(potts_folder_2, DUMMY_SEQUENCE_FILE_2, "one_submat")
        output_folder = pathlib.Path(tempfile.mkdtemp())
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_objects([obj1,obj2], output_folder)
        output_file = '/tmp/'+next(tempfile._get_candidate_names())
        get_seqs_aligned_in_fasta_file(aligned_positions_with_gaps, [obj1,obj2], output_file)
        records = list(SeqIO.parse(str(output_file), "fasta"))
        assert(str(records[0].seq)=="C-WP")
        assert(str(records[1].seq)=="CGWP")



if __name__=='__main__':
    unittest.main()



