import unittest
import shutil, tempfile
import pathlib
import numpy as np

from makepotts.potts_model import *
from ppalign.compute_scores import *
from ppalign.call_ppalign import *

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
        edges_map = get_edges_map(self.mrf, 100)
        selfcomp = compute_selfscore(self.mrf, edges_map, use_v=True, use_w=False)
        aln_dict, infos_solver = align_two_potts_models([self.mrf, self.mrf], self.output_folder, use_w=False)
        UB = infos_solver['UB']
        self.assertTrue(abs(selfcomp-UB)<1)


    def test_selfscore_no_v(self):
        edges_map = get_edges_map(self.mrf, 100)
        selfcomp = compute_selfscore(self.mrf, edges_map, use_v=False, use_w=True)
        aln_dict, infos_solver = align_two_potts_models([self.mrf, self.mrf], self.output_folder, use_v=False)
        UB = infos_solver['UB']
        self.assertTrue(abs(selfcomp-UB)<1)


    def test_selfscore(self):
        v_score_function = scalar_product
        w_score_function = scalar_product
        edges_map = get_edges_map(self.mrf, 100)
        selfcomp = compute_selfscore(self.mrf, edges_map, use_v=True, use_w=True)
        aln_dict, infos_solver = align_two_potts_models([self.mrf, self.mrf], self.output_folder, use_w=True)
        UB = infos_solver['UB']
        self.assertTrue(abs(selfcomp-UB)<1)


    def test_selfscore_alpha_w(self):
        edges_map = get_edges_map(self.mrf, 100)
        selfcomp_normal = compute_selfscore(self.mrf, edges_map, use_v=False, use_w=True)
        alpha_w = 2
        selfcomp_alpha = compute_selfscore(self.mrf, edges_map, use_v=False, use_w=True, alpha_w=alpha_w)
        assert(selfcomp_alpha==alpha_w*selfcomp_normal)

    def test_background_v0(self):
        v0 = get_background_v0("identity")
        assert(sum(v0)<0.00001)


    def test_get_gap_cost(self):
        gap_open=14
        assert(get_total_gap_cost({"pos_ref":[4,5,7,8,9,10], "pos_2":[0,1,2,3,4,5]}, gap_open)==gap_open)
        assert(get_total_gap_cost({"pos_ref":[4,5,7,8,10], "pos_2":[0,1,2,3,5]}, gap_open)==2*gap_open)
        assert(get_total_gap_cost({"pos_ref":[4,5,7,8,9], "pos_2":[0,1,2,3,4]}, gap_open)==gap_open)


    def test_get_total_score(self):
        params = {'alpha_w':1, 'remove_v0':True, 'offset_v':1, 'gap_open':10, 'v_rescaling_function':'identity'}
        aln_dict, infos_solver = align_two_potts_models([self.mrf, self.mrf], self.output_folder, use_w=False, **params)
        LB = infos_solver['LB']
        total_score = get_score_for_alignment([self.mrf, self.mrf], aln_dict, **params)
        assert((LB-total_score)<0.1)




if __name__=='__main__':
    unittest.main()
