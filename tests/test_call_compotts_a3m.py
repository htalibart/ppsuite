import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

from compotts.compotts_object import *
from compotts.call_compotts import *

import pkg_resources
EXAMPLES_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/test_call_compotts_hhblits_object/'))



class Test_Call_ComPotts_a3m(unittest.TestCase):

    def setUp(self):
        PROTEIN_NAME = "1cc8"
        seq_file = EXAMPLES_FOLDER/(PROTEIN_NAME+".fasta")
        a3m_file = EXAMPLES_FOLDER/(PROTEIN_NAME+".a3m")
        self.output_folder = pathlib.Path(tempfile.mkdtemp())
        self.obj = ComPotts_Object(seq_file=seq_file, a3m_file=a3m_file, input_folder=self.output_folder)

    def tearDown(self):
        shutil.rmtree(self.output_folder)


    def test_align_compotts_object_to_itself(self):
        aligned_positions, infos_solver = align_two_objects([self.obj, self.obj], self.output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)


if __name__=='__main__':
    unittest.main()
