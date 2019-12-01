import unittest
import shutil, tempfile
import pathlib
import numpy as np

from comutils.potts_model import *
from compotts.compute_scores import *
from compotts.call_compotts import *

import pkg_resources
from resources_manager import *

class Test_Compute_Scores(unittest.TestCase):

    def setUp(self):
        self.mrf = Potts_Model.from_msgpack(MRF_1CC8)
        self.output_folder = pathlib.Path(tempfile.mkdtemp())


    def tearDown(self):
        shutil.rmtree(self.output_folder)


    def test_selfscore_no_w(self):
        v_score_function = scalar_product
        w_score_function = scalar_product
        edges_map = get_edges_map(self.mrf, "")
        selfcomp = compute_selfscore(self.mrf, edges_map, v_score_function, w_score_function, use_v=True, use_w=False)
        aln_dict, infos_solver = align_two_potts_models([self.mrf, self.mrf], self.output_folder, use_w=False)
        UB = infos_solver['UB']
        self.assertTrue(abs(selfcomp-UB)<1)


    def test_selfscore(self):
        v_score_function = scalar_product
        w_score_function = scalar_product
        edges_map = get_edges_map(self.mrf, "")
        selfcomp = compute_selfscore(self.mrf, edges_map, v_score_function, w_score_function, use_v=True, use_w=True)
        aln_dict, infos_solver = align_two_potts_models([self.mrf, self.mrf], self.output_folder, use_w=True)
        UB = infos_solver['UB']
        self.assertTrue(abs(selfcomp-UB)<1)


if __name__=='__main__':
    unittest.main()
