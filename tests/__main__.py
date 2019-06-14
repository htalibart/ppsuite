import unittest
from compotts_wrapper.compotts_object import *
from compotts_wrapper.call_compotts import *


TEST_OUTPUT_FOLDER = "tests_output/"
EXAMPLES_FOLDER = "examples/"
SIMPLE_TEST = "1cc8"

class TestComPottsObject(unittest.TestCase):

    def test_dummy(self):
        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")

    def test_import_hhblits(self):
        test_name = EXAMPLES_FOLDER+SIMPLE_TEST
        a3m_file = test_name+".a3m"
        seq_file = test_name+".fasta"
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, output_folder=TEST_OUTPUT_FOLDER, nb_sequences=200)

    def test_to_one_hot(self):
        seq_file = EXAMPLES_FOLDER+SIMPLE_TEST+".fasta" 
        obj = ComPotts_Object.from_seq_file_to_one_hot(seq_file, output_folder=TEST_OUTPUT_FOLDER)


    def test_align_compotts_object_to_itself(self):
        test_name = EXAMPLES_FOLDER+SIMPLE_TEST
        a3m_file = test_name+".a3m"
        seq_file = test_name+".fasta"
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, output_folder=TEST_OUTPUT_FOLDER, nb_sequences=200)
        aln_res_file = TEST_OUTPUT_FOLDER+SIMPLE_TEST+"_"+SIMPLE_TEST+"_aln.csv"
        info_res_file = TEST_OUTPUT_FOLDER+SIMPLE_TEST+"_"+SIMPLE_TEST+"_info.csv"
        align_two_objects([obj, obj], aln_res_file, info_res_file, scores_folder=TEST_OUTPUT_FOLDER)

        

if __name__=='__main__':
    unittest.main()
