import unittest

import numpy as np

from compotts_wrapper.compotts_object import *
from compotts_wrapper.call_compotts import *
from compotts_wrapper.manage_positions import *
from compotts_wrapper.rescaling import *
import create_fake_data as crfake
from potts_model import *
import files_management as fm

TEST_OUTPUT_FOLDER = "tests_output/"
EXAMPLES_FOLDER = "examples/"
SIMPLE_TEST = "1cc8"

def are_templates_aligned(template2, aligned_positions):
    good_alignment = True
    pos2 = 0
    while (good_alignment) and (pos2<len(template2)):
        t = template2[pos2]
        if t[0]=='[':
            pos1 = int(t[1:len(t)-1].replace('-',''))
            real_pos_2 = get_pos_aligned_at_pos(aligned_positions, pos1)
            good_alignment = (pos2==real_pos_2)
        pos2+=1
    return good_alignment



class TestComPottsObject(unittest.TestCase):

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
        output_folder = fm.create_folder(TEST_OUTPUT_FOLDER+SIMPLE_TEST+"_"+SIMPLE_TEST+"/")
        a3m_file = test_name+".a3m"
        seq_file = test_name+".fasta"
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, output_folder=output_folder, nb_sequences=200)
        aligned_positions, infos_solver = align_two_objects([obj, obj], output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)


    def test_align_small_fake_mrfs(self):
        output_folder = fm.create_folder(TEST_OUTPUT_FOLDER+"fake/")
        templates = [["1", "0", "3", "2"],["[0]", "y", "[1]", "[2]"]]
        alnfnames = [output_folder+"fake_"+str(i)+".aln" for i in range(2)]
        fastafnames = [output_folder+"fake_"+str(i)+".fasta" for i in range(2)]
        crfake.main(templates, alnfnames, fastafnames)
        mrfs = [Potts_Model.from_training_set(fastafnames[i], output_folder+"fake_"+str(i)+".mrf") for i in range(2)]
        aligned_positions, infos_solver = align_two_potts_models(mrfs, output_folder)
        self.assertTrue(are_templates_aligned(templates[1], aligned_positions))


    def test_align_one_hot_to_itself(self):
        test_name = SIMPLE_TEST+"_one_hot"
        output_folder = fm.create_folder(TEST_OUTPUT_FOLDER+test_name+"/")
        seq_file = EXAMPLES_FOLDER+SIMPLE_TEST+".fasta"
        aligned_positions, infos_solver = align_one_hot([seq_file, seq_file], output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)


    def test_identity_rescaling(self):
        for x in [12, [1, 2, 3], np.ones((3,3))]:
            self.assertTrue(np.array_equal(identity(x), x))


    def test_count_edges(self):
        self.assertEqual(count_edges(np.ones((3,3))),9)
        self.assertEqual(count_edges(np.zeros((3,3))),0)


    def test_get_edges_map(self):
        mrf = Potts_Model.from_msgpack(EXAMPLES_FOLDER+SIMPLE_TEST+".mrf")
        edges_map = get_edges_map(mrf, 0)
        self.assertEqual(sum([abs(edges_map[i][i]) for i in range(mrf.ncol)]), 0)


    def test_get_seqs_aligned(self):
        objs = []
        for k in range(2):
            seq_file = EXAMPLES_FOLDER+"fake_seq_"+str(k)+".fasta"
            objs.append(ComPotts_Object.from_seq_file_to_one_hot(seq_file, output_folder=TEST_OUTPUT_FOLDER))
        aligned_positions = {"pos_ref":[4,5,7,8,9,10], "pos_2":[0,1,2,3,4,5]} 
        seqs_aligned = get_seqs_aligned(aligned_positions, objs)
        self.assertEqual(seqs_aligned[0], '----MA-IKEH-')
        self.assertEqual(seqs_aligned[1], '----MA-IKDH-')


    def test_rescale_mrf(self):
        mrf = Potts_Model.from_msgpack(EXAMPLES_FOLDER+SIMPLE_TEST+".mrf")
        resc_mrf = get_rescaled_mrf(mrf, "original_rescaling")


    def test_rescale_compotts_object(self):
        test_name = EXAMPLES_FOLDER+SIMPLE_TEST
        a3m_file = test_name+".a3m"
        seq_file = test_name+".fasta"
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, output_folder=TEST_OUTPUT_FOLDER, nb_sequences=200, rescaling_function="original_rescaling")


    def test_aligned_rescaled_mrf_to_itself(self):
        test_name = EXAMPLES_FOLDER+SIMPLE_TEST
        rescaling_function="original_rescaling"
        output_folder = fm.create_folder(TEST_OUTPUT_FOLDER+SIMPLE_TEST+"_"+SIMPLE_TEST+"_"+rescaling_function+"/")
        a3m_file = test_name+".a3m"
        seq_file = test_name+".fasta"
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, output_folder=output_folder, nb_sequences=200, rescaling_function=rescaling_function)
        aligned_positions, infos_solver = align_two_objects([obj, obj], output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)



if __name__=='__main__':
    unittest.main()
