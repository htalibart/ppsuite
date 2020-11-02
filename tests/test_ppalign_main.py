import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

import pkg_resources
from tests.resources_manager import *

from ppalign.__main__ import *


class Test_PPalign_Main(unittest.TestCase):

    def setUp(self):
        self.output_folder = pathlib.Path(tempfile.mkdtemp())
        self.feature_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        self.feature_folder_1 = self.feature_folder
        self.feature_folder_2 = self.feature_folder
        shutil.copytree(FEATURE_FOLDER, self.feature_folder)

    def tearDown(self):
        shutil.rmtree(self.output_folder)
        shutil.rmtree(self.feature_folder)

    def test_self_alignment_from_folders(self):
        ppalign_args = ["--feature_folder_1", str(self.feature_folder_1), "--feature_folder_2", str(self.feature_folder_1), "--output_folder", str(self.output_folder), "--get_training_sets_fasta_aln", "--get_sequences_fasta_aln"]
        res_ppalign = main(ppalign_args)
        assert(res_ppalign['infos_solver']['UB']==res_ppalign['infos_solver']['selfcomp1'])

    def test_self_alignment_from_potts_model_files(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder)]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)
        assert(res_ppalign['infos_solver']['UB']==res_ppalign['infos_solver']['selfcomp1'])

    def test_exponential_scalar_product(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--exp"]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)


    def test_offset_v(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--offset", str(2)]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)

    def test_remove_v0(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--remove_v0"]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)

    def test_remove_v0_rescaled(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--remove_v0", "--rescale_removed_v0", "--v_rescaling_function", "simulate_uniform_pc_on_v"]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)

    def test_remove_negative_couplings(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--w_rescaling_function", "remove_negative_couplings"]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)


    def test_perc_w(self):
        ppalign_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--w_percent", str(10)]
        res_ppalign = main(ppalign_args)
        print(res_ppalign)



if __name__=='__main__':
    unittest.main()
