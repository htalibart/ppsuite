import unittest
import shutil, tempfile
import pathlib
import os
import numpy as np

from compotts.compotts_object import *
from compotts.call_compotts import *


import pkg_resources
EXAMPLES_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/test_call_compotts_one_hot/'))


class Test_Call_ComPotts_OneHot(unittest.TestCase):

    def setUp(self):
        PROTEIN_NAME = "1cc8"
        seq_file = EXAMPLES_FOLDER/(PROTEIN_NAME+".fasta")
        self.output_folder = pathlib.Path(tempfile.mkdtemp())
        self.obj = ComPotts_Object(sequence_file=seq_file, input_folder=self.output_folder, mrf_type="one_hot")

    def tearDown(self):
        shutil.rmtree(self.output_folder)

    def test_align_one_hot_to_itself(self):
        aligned_positions, infos_solver = align_two_objects([self.obj, self.obj], self.output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)


if __name__=='__main__':
    unittest.main()
