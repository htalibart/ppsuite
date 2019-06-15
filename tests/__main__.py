import unittest
from compotts_wrapper.compotts_object import *
from compotts_wrapper.call_compotts import *
from compotts_wrapper.manage_positions import *
from compotts_wrapper.get_infos import *
import create_fake_data as crfake
from potts_model import *

TEST_OUTPUT_FOLDER = "tests_output/"
EXAMPLES_FOLDER = "examples/"
SIMPLE_TEST = "1cc8"

def are_templates_aligned(template2, aln_res_file):
    good_alignment = True
    pos2 = 0
    while (good_alignment) and (pos2<len(template2)):
        t = template2[pos2]
        if t[0]=='[':
            pos1 = int(t[1:len(t)-1].replace('-',''))
            real_pos_2 = get_pos_aligned_at_pos(aln_res_file, pos1)
            good_alignment = (pos2==real_pos_2)
        pos2+=1
    return good_alignment


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
        similarity_global = get_info("similarity_global", info_res_file)
        assert(similarity_global==1)



    def test_align_small_fake_mrfs(self):
        f = TEST_OUTPUT_FOLDER
        templates = [["1", "0", "3", "2"],["[0]", "y", "[1]", "[2]"]]
        alnfnames = [f+"fake_"+str(i)+".aln" for i in range(2)]
        fastafnames = [f+"fake_"+str(i)+".fasta" for i in range(2)]
        crfake.main(templates, alnfnames, fastafnames)
        mrfs = [Potts_Model.from_training_set(fastafnames[i], f+"fake_"+str(i)+".mrf") for i in range(2)]
        aln_res_file = f+"fake_aln.csv"
        info_res_file = f+"fake_info.csv"
        align_two_potts_models(mrfs, aln_res_file, info_res_file, output_folder=TEST_OUTPUT_FOLDER)
        assert(are_templates_aligned(templates[1], aln_res_file))

if __name__=='__main__':
    unittest.main()
