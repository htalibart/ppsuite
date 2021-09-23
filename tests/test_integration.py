import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

import pkg_resources
from tests.resources_manager import *

from ppalign.__main__ import main as ppalign_main
from makepotts.potts_object import main as potts_main


class Test_Integration(unittest.TestCase):

    def setUp(self):
        self.output_folder = pathlib.Path(tempfile.mkdtemp())
        self.potts_folder_1 = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        self.potts_folder_2 = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        self.potts_folder_3 = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))

    def tearDown(self):
        pass
        for folder in [self.output_folder, self.potts_folder_1, self.potts_folder_2, self.potts_folder_3]:
            if folder.is_dir():
                shutil.rmtree(folder)

    def test_insertion(self):
        potts_main(["--potts_folder", str(self.potts_folder_1), "--aln_with_insertions", str(INSERTION_RESOURCES_FOLDER/'ACCCD.a3m'), "--use_insertion_penalties", "--dont_trim_alignment"])
        potts_main(["--potts_folder", str(self.potts_folder_1), "--aln_with_insertions", str(INSERTION_RESOURCES_FOLDER/'ACCCCD.a3m'), "--use_insertion_penalties", "--dont_trim_alignment", "--dont_filter_alignment"])
        potts_main(["--potts_folder", str(self.potts_folder_2), "--aln_with_insertions", str(INSERTION_RESOURCES_FOLDER/'AiD.a3m'), "--use_insertion_penalties", "--dont_trim_alignment", "--dont_filter_alignment"])
        potts_main(["--potts_folder", str(self.potts_folder_3), "--aln_with_insertions", str(INSERTION_RESOURCES_FOLDER/'AD.a3m'), "--use_insertion_penalties", "--dont_trim_alignment"])
        res_ppalign_1 = ppalign_main(["--potts_folder_1", str(self.potts_folder_1), "--potts_folder_2", str(self.potts_folder_2), "--output_folder", str(self.output_folder), "--get_sequences_fasta_aln", "--use_insertion_penalties", "--t_limit", str(1)])
        res_ppalign_2 = ppalign_main(["--potts_folder_1", str(self.potts_folder_1), "--potts_folder_2", str(self.potts_folder_3), "--output_folder", str(self.output_folder), "--get_sequences_fasta_aln", "--use_insertion_penalties"])
        assert(res_ppalign_1['infos_solver']['UB']>res_ppalign_2['infos_solver']['UB'])


    def test_insertion_coeff(self):
        potts_main(["--potts_folder", str(self.potts_folder_1), "--aln_with_insertions", str(INSERTION_RESOURCES_FOLDER/'ACCCD.a3m'), "--use_insertion_penalties", "--dont_trim_alignment"])
        potts_main(["--potts_folder", str(self.potts_folder_2), "--aln_with_insertions", str(INSERTION_RESOURCES_FOLDER/'AiD.a3m'), "--use_insertion_penalties", "--dont_trim_alignment", "--dont_filter_alignment"])
        res_ppalign_1 = ppalign_main(["--potts_folder_1", str(self.potts_folder_1), "--potts_folder_2", str(self.potts_folder_2), "--output_folder", str(self.output_folder), "--get_sequences_fasta_aln", "--use_insertion_penalties", "--t_limit", str(1)])
        res_ppalign_2 = ppalign_main(["--potts_folder_1", str(self.potts_folder_1), "--potts_folder_2", str(self.potts_folder_2), "--output_folder", str(self.output_folder), "--get_sequences_fasta_aln", "--use_insertion_penalties", "--t_limit", str(1), "--insertion_penalties_coefficient", str(10)])

        assert(res_ppalign_1['infos_solver']['LB']>res_ppalign_2['infos_solver']['LB'])




if __name__=='__main__':
    unittest.main()
