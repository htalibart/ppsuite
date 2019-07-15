import unittest
import shutil, tempfile

import numpy as np

from compotts.compotts_object import *
from compotts.manage_positions import *


import pkg_resources
EXAMPLES_FOLDER = pkg_resources.resource_filename(__name__,'examples/test_manage_positions/')

class Test_ManagePositions(unittest.TestCase):

    def setUp(self):
        self.output_folder = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.output_folder)

    def test_get_seqs_aligned(self):
        objs = []
        for k in range(2):
            seq_file = EXAMPLES_FOLDER+"fake_seq_"+str(k)+".fasta"
            objs.append(ComPotts_Object.from_seq_file_to_one_hot(seq_file, self.output_folder))
        aligned_positions = {"pos_ref":[4,5,7,8,9,10], "pos_2":[0,1,2,3,4,5]}
        seqs_aligned = get_seqs_aligned(aligned_positions, objs)
        self.assertEqual(seqs_aligned[0], '----MA-IKEH-')
        self.assertEqual(seqs_aligned[1], '----MA-IKDH-')



if __name__=='__main__':
    unittest.main()
