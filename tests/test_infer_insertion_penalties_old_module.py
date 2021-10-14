import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils import files_management as fm
from makepotts.handle_insertions import *
from makepotts.potts_object import *
from makepotts.handle_insertions_julia import *


import time

class Test_Infer_Insertion_Penalties(unittest.TestCase):

    def setUp(self):
        self.potts_folder = pathlib.Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.potts_folder)
        #pass

#    def test_count_insertions(self):
#       a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"count_insertions_example.a3m")
#       ins = count_insertions(a3m_file)
#       assert(np.array_equiv(ins,np.array([[0,0,0,0],[0,1,0,0],[0,3,0,0]])))
#
#    def test_infer_insertion_penalties(self):
#        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"count_insertions_example.a3m")
#        insertion_penalties = infer_insertion_penalties(a3m_file)
#
#    def test_lower_trimmed_columns_insertions(self):
#        a3m_file = INSERTION_RESOURCES_FOLDER/'Ac-D.a3m'
#        output_file = '/tmp/test_a3m_trim.a3m'
#        lower_case_trimmed_columns(a3m_file, output_file, [0,2])
#        records = list(SeqIO.parse(str(output_file), 'fasta'))
#        assert(str(records[2].seq)=='AcfD')
#        os.remove(str(output_file))
#
#    def test_lower_trimmed_columns_insertions_full_object(self):
#        a3m_file = INSERTION_RESOURCES_FOLDER/'Ac-D.a3m'
#        po = Potts_Object.from_hhblits_files(self.potts_folder, aln_with_insertions=a3m_file, use_insertion_penalties=True, filter_alignment=False, trim_alignment=True, trimal_gt=0.9)
#        assert(po.potts_model.ncol==2)
#        assert(po.insertion_penalties['open'][1]<po.insertion_penalties['open'][2])
# 
#    def test_infer_insertion_penalties_with_pseudocounts(self):
#        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"count_insertions_example.a3m")
#        insertion_penalties = infer_insertion_penalties(a3m_file, pc_insertions_tau=0)
#        insertion_penalties_with_pseudocounts = infer_insertion_penalties(a3m_file, pc_insertions_tau=0.1)
#        print(insertion_penalties)
#        print(insertion_penalties_with_pseudocounts)
#        assert(insertion_penalties["open"][1]<insertion_penalties_with_pseudocounts["open"][1])
#        assert(insertion_penalties["extend"][1]<insertion_penalties_with_pseudocounts["extend"][1])
#
#    def test_infer_insertion_penalties_not_negative(self):
#        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"gap_open.a3m")
#        insertion_penalties_without_pc = infer_insertion_penalties(a3m_file, pc_insertions_tau=0)
#        print(insertion_penalties_without_pc)
#        assert(insertion_penalties_without_pc["open"][1]<0)
#        insertion_penalties_with_pc = infer_insertion_penalties(a3m_file, pc_insertions_tau=0.5)
#        print(insertion_penalties_with_pc)
#        assert(insertion_penalties_with_pc["open"][1]>0)


#    def test_learning_coeff_insertions(self):
#        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER)/"Q9AZ42_filtered.a3m"
#        insertion_penalties_new_learning_rate = infer_insertion_penalties(a3m_file, pc_insertions_tau=0, learning_coeff_insertions=1)
#        print(insertion_penalties_new_learning_rate)
#        insertion_penalties_standard_learning_rate = infer_insertion_penalties(a3m_file, pc_insertions_tau=0, learning_coeff_insertions=1e-3)
#        print(insertion_penalties_standard_learning_rate)
#        L = get_length_ins_file(a3m_file) 
#        for col in range(L+1):
#            for gap_type in ['open','extend']:
#                assert(insertion_penalties_standard_learning_rate[gap_type][col]-insertion_penalties_new_learning_rate[gap_type][col]<0.001)
#


    def test_compare_with_julia(self):
        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER)/"Q9AZ42_filtered.a3m"
        seed_length = get_length_ins_file(a3m_file)
        output_julia_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        start = time.time()
        infer_insertion_penalties_in_file_using_dcabuild(a3m_file, seed_length, output_julia_file)
        end_julia = time.time()
        time_julia = end_julia-start
        insertion_penalties_julia = get_insertion_penalties_from_file(output_julia_file)
        start = time.time()
        #insertion_penalties = infer_insertion_penalties(a3m_file, maxit_infer_insertions=6e7, tol_infer_insertions=1e-3, learning_coeff_insertions=1e-3, freq_insert_min=1e-3, pc_insertions_tau=0)
        insertion_penalties = infer_insertion_penalties(a3m_file)
        end = time.time()
        time_pp = end-start
        print("time_julia=",time_julia, "time_pp=",time_pp)
        print(insertion_penalties_julia)
        print(insertion_penalties)
        output_julia_file.unlink()



if __name__=='__main__':
    unittest.main()
