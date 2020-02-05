import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

import pkg_resources
from tests.resources_manager import *

from compotts.__main__ import *


class Test_Call_ComPotts(unittest.TestCase):

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
        compotts_args = ["--feature_folder_1", str(self.feature_folder_1), "--feature_folder_2", str(self.feature_folder_1), "--output_folder", str(self.output_folder), "--get_training_sets_fasta_aln", "--get_sequences_fasta_aln"]
        res_compotts = main(compotts_args)
        assert(res_compotts['infos_solver']['UB']==res_compotts['infos_solver']['selfcomp1'])

    def test_self_alignment_from_potts_model_files(self):
        compotts_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder)]
        res_compotts = main(compotts_args)
        print(res_compotts)
        assert(res_compotts['infos_solver']['UB']==res_compotts['infos_solver']['selfcomp1'])

    def test_self_alignment_with_rescaling(self):
        compotts_args = ["--potts_model_file_1", str(self.feature_folder_1/"potts_model.mrf"), "--potts_model_file_2", str(self.feature_folder_1/"potts_model.mrf"), "--output_folder", str(self.output_folder), "--rescaling_function", "add_number", "--shift", str(3)]
        res_compotts = main(compotts_args)
        print(res_compotts)


if __name__=='__main__':
    unittest.main()
