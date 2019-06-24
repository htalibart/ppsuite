import unittest
import shutil, tempfile
import os
import numpy as np

from compotts.compotts_object import *
from compotts.call_compotts import *


import pkg_resources
EXAMPLES_FOLDER = pkg_resources.resource_filename(__name__,'examples/test_call_compotts_one_hot/')


class Test_Call_ComPotts_OneSeq(unittest.TestCase):

    def setUp(self):
        PROTEIN_NAME = "1cc8"
        seq_file = os.path.join(EXAMPLES_FOLDER,PROTEIN_NAME+".fasta")
        self.output_folder = tempfile.mkdtemp()
        self.obj = ComPotts_Object.from_seq_file_via_ccmpred(seq_file, self.output_folder)

    def tearDown(self):
        shutil.rmtree(self.output_folder)

    def test_align_seq_to_itself(self):
        aligned_positions, infos_solver = align_two_objects([self.obj, self.obj], self.output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)


if __name__=='__main__':
    unittest.main()
