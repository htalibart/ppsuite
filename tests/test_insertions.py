import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

from makepotts.potts_object import *
from ppalign.call_ppalign import *

import pkg_resources
from tests.resources_manager import *

import comutils.create_fake_data as crfake

from vizpotts import vizpotts #TEMP


def are_templates_aligned(template2, aligned_positions):
    good_alignment = True
    pos2 = 0
    while (good_alignment) and (pos2<len(template2)):
        t = template2[pos2]
        if t[0]=='[':
            pos1 = int(t[1:len(t)-1].replace('-',''))
            real_pos_2 = get_pos_aligned_at_pos(aligned_positions, pos1)
            good_alignment = (pos2==real_pos_2)
        pos2+=1
    return good_alignment

def get_vi_conserved_letter(conserved_letter, via=1):
    via_conserved = via
    vi = np.ones(21)*(-via_conserved/(21-1))
    vi[conserved_letter] = via_conserved
    return vi


def get_fake_model(conserved_letters_list, ijabs=[], wijab=10, via=1):
    v = np.array([get_vi_conserved_letter(c, via=via) for c in conserved_letters_list])
    w = np.zeros((len(conserved_letters_list), len(conserved_letters_list), 21, 21))
    for ijab in ijabs:
        w[ijab] = wijab
        w[ijab[1],ijab[0],ijab[2],ijab[3]] = wijab
    return Potts_Model.from_parameters(v,w)


class Test_Insertions(unittest.TestCase):

    def setUp(self):
        self.output_folder = pathlib.Path(tempfile.mkdtemp())
        self.potts_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.potts_folder)
        self.object = Potts_Object.from_folder(potts_folder=self.potts_folder)
        self.potts_model = self.object.potts_model

    def tearDown(self):
        shutil.rmtree(self.output_folder)
        shutil.rmtree(self.potts_folder)


    def test_costly_open(self):
        mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
        insert_costs = [{"open":[0,1000,0,0],"extend":[0,0,0,0]}, {"open":[0,0,0,0,0,0],"extend":[0,0,0,0,0,0]}]
        expected_aligned_positions = {"pos_ref":[1,2], "pos_2":[3,4]}
        expected_aligned_positions_with_gaps = {"pos_ref":['-','-','-',0,1,2], "pos_2":[0,1,2,'-',3,4]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100, t_limit=100)
        self.assertEqual(aligned_positions, expected_aligned_positions)
        self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_not_costly_open(self):
        mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
        insert_costs = [{"open":[0,0,0,0],"extend":[0,0,0,0]}, {"open":[0,0,0,0,0,0],"extend":[0,0,0,0,0,0]}]
        expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,3,4]}
        expected_aligned_positions_with_gaps = {"pos_ref":[0,'-','-',1,2], "pos_2":[0,1,2,3,4]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100, t_limit=100)
        self.assertEqual(aligned_positions, expected_aligned_positions)
        self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])



    def test_costly_open_other_model(self):
        mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
        insert_costs = [{"open":[0,0,0,0],"extend":[0,0,0,0]}, {"open":[1000,1000,1000,1000,1000,1000],"extend":[0,0,0,0,0,0]}]
        expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,3,4]}
        expected_aligned_positions_with_gaps = {"pos_ref":[0,'-','-',1,2], "pos_2":[0,1,2,3,4]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(aligned_positions, expected_aligned_positions)
        self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_not_very_costly_extend_1(self):
       mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
       insert_costs = [{"open":[0,0,0,0],"extend":[0,0.33333,0,0]}, {"open":[0,0,0,0,0,0],"extend":[0,0,0,0,0,0]}]
       expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,3,4]}
       expected_aligned_positions_with_gaps = {"pos_ref":[0,'-','-',1,2], "pos_2":[0,1,2,3,4]}
       aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100)
       self.assertEqual(aligned_positions, expected_aligned_positions)
       self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
       self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_not_very_costly_extend_longer(self):
       mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,5,1,2])]
       insert_costs = [{"open":[0,0,0,0],"extend":[0,0.33333,0,0]}, {"open":[0,0,0,0,0,0,0],"extend":[0,0,0,0,0,0,0]}]
       expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,4,5]}
       expected_aligned_positions_with_gaps = {"pos_ref":[0,'-','-','-',1,2], "pos_2":[0,1,2,3,4,5]}
       aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100)
       self.assertEqual(aligned_positions, expected_aligned_positions)
       self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
       self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_not_very_costly_extend_with_open(self):
       mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
       insert_costs = [{"open":[0,0.1,0,0],"extend":[0,0.33333,0,0]}, {"open":[0,0,0,0,0,0],"extend":[0,0,0,0,0,0]}]
       expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,3,4]}
       expected_aligned_positions_with_gaps = {"pos_ref":[0,'-','-',1,2], "pos_2":[0,1,2,3,4]}
       aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100)
       self.assertEqual(aligned_positions, expected_aligned_positions)
       self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
       self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_not_very_costly_extend_and_useless(self):
       mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,1,2])]
       insert_costs = [{"open":[0,0,0,0],"extend":[0,0.33333,0,0]}, {"open":[0,0,0,0,0],"extend":[0,0,0,0,0]}]
       expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,2,3]}
       expected_aligned_positions_with_gaps = {"pos_ref":[0,'-',1,2], "pos_2":[0,1,2,3]}
       aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100)
       self.assertEqual(aligned_positions, expected_aligned_positions)
       self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
       self.assertEqual(infos_solver['UB'], infos_solver['LB'])



    def test_costly_extend_1(self):
       mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
       insert_costs = [{"open":[0,0,0,0],"extend":[0,1000,0,0]}, {"open":[0,0,0,0,0,0],"extend":[0,0,0,0,0,0]}]
       expected_aligned_positions = {"pos_ref":[1,2], "pos_2":[3,4]}
       aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=-100)
       self.assertEqual(aligned_positions, expected_aligned_positions)
       self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_costly_extend_other_model(self):
       mrfs = [get_fake_model([0,1,2]), get_fake_model([0,3,4,1,2])]
       insert_costs = [{"open":[0,0,0,0],"extend":[0,0,0,0]}, {"open":[0,0,0,0,0,0],"extend":[1000,1000,1000,1000,1000,1000]}]
       expected_aligned_positions = {"pos_ref":[0,1,2], "pos_2":[0,3,4]}
       expected_aligned_positions_with_gaps = {"pos_ref":[0,'-','-',1,2], "pos_2":[0,1,2,3,4]}
       aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
       self.assertEqual(aligned_positions, expected_aligned_positions)
       self.assertEqual(aligned_positions_with_gaps, expected_aligned_positions_with_gaps)
       self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_costly_extend_2(self):
        mrfs = [get_fake_model([0,1,2,3]), get_fake_model([0,4,5,6,1,2,3])]
        insert_costs = [{"open":[0,0,0,0,0],"extend":[0,1000,0,0,0]}, {"open":[0,0,0,0,0,0,0,0],"extend":[0,0,0,0,0,0,0,0]}]
        expected_aligned_positions = {"pos_ref":[1,2,3], "pos_2":[4,5,6]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(aligned_positions, expected_aligned_positions)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])

    def test_no_external_cost(self):
        nb_before = 1
        nb_after = 1
        mrfs = [get_fake_model([0]*nb_before+[1,2]+[0]*nb_after), get_fake_model([1,2])]
        insert_costs = [{"open":[0]*(mrfs[0].ncol+1), "extend":[0]+[0]*(mrfs[0].ncol+1)},
                        {"open":[0,1000,0], "extend":[0,0,0]}
                        ]
        expected_aligned_positions = {"pos_ref":[nb_before,nb_before+1], "pos_2":[0,1]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(aligned_positions, expected_aligned_positions)

    def test_external_cost(self):
        nb_before = 1
        nb_after = 1
        mrfs = [get_fake_model([0]*nb_before+[1,2]+[0]*nb_after), get_fake_model([1,2])]
        insert_costs = [{"open":[0]*(mrfs[0].ncol+1), "extend":[0]+[0]*(mrfs[0].ncol+1)},
                        {"open":[1000,0,300], "extend":[0,0,0]}
                        ]
        not_expected_aligned_positions = {"pos_ref":[nb_before,nb_before+1], "pos_2":[0,1]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001, n_limit_param=1)
        self.assertNotEqual(aligned_positions, not_expected_aligned_positions)


    def test_external_cost_after(self):
        open_after = 0.33333
        mrfs = [get_fake_model([0,1]), get_fake_model([0])]
        insert_costs = [{"open":[0]*(mrfs[0].ncol+1), "extend":[0]*(mrfs[0].ncol+1)},
                        {"open":[0,open_after], "extend":[0]*(mrfs[1].ncol+1)}
                        ]
        expected_aligned_positions = {"pos_ref":[0], "pos_2":[0]}
        expected_LB = scalar_product(mrfs[0].v[0],mrfs[1].v[0])-open_after
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001, n_limit_param=1)
        self.assertEqual(aligned_positions, expected_aligned_positions)
        assert(abs(infos_solver['LB']-expected_LB)<=0.001)


    def test_external_cost_before(self):
        open_before = 0.33333
        mrfs = [get_fake_model([1,0]), get_fake_model([0])]
        insert_costs = [{"open":[0]*(mrfs[0].ncol+1), "extend":[0]*(mrfs[0].ncol+1)},
                        {"open":[open_before,0], "extend":[0]*(mrfs[1].ncol+1)}
                        ]
        expected_aligned_positions = {"pos_ref":[1], "pos_2":[0]}
        expected_LB = scalar_product(mrfs[1].v[0],mrfs[1].v[0])-open_before
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001, n_limit_param=1)
        self.assertEqual(aligned_positions, expected_aligned_positions)
        print(expected_LB, infos_solver['LB'])
        assert(abs(infos_solver['LB']-expected_LB)<=0.001)


#    def test_w_score_open_gap(self): # WILL NOT CONVERGE
#        gap_open=1000
#        gap_extend=0
#        wijab=0.25
#        via=1
#        mrf1 = get_fake_model([0,1,2], ijabs=[(0,2,0,0)], wijab=wijab, via=via)
#        mrf2 = get_fake_model([0,2], ijabs=[(0,1,0,0)], wijab=wijab, via=via)
#        insert_costs = [{"open":[0]*(mrf1.ncol+1), "extend":[0]*(mrf1.ncol+1)},
#                        {"open":[0,gap_open,0], "extend":[0,gap_extend,0]}
#                        ]
#        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models([mrf1,mrf2], self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001, n_limit_param=1, iter_limit_param=1000000000, theta=2)
#        self.assertEqual(infos_solver['UB'], infos_solver['LB'])
#        not_expected_aligned_positions = {"pos_ref":[0,2], "pos_2":[0,1]}
#        self.assertNotEqual(aligned_positions, not_expected_aligned_positions)

    def test_w_score_open_gap_without_w(self):
        gap_open=1000
        gap_extend=0
        mrf1 = get_fake_model([0,1,2])
        mrf2 = get_fake_model([0,2])
        insert_costs = [{"open":[0]*(mrf1.ncol+1), "extend":[0]*(mrf1.ncol+1)},
                        {"open":[0,gap_open,0], "extend":[0,gap_extend,0]}
                        ]
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models([mrf1,mrf2], self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])
        possible_expected_aligned_positions = [{"pos_ref":[0], "pos_2":[0]}, {"pos_ref":[2], "pos_2":[1]}]
        assert((aligned_positions==possible_expected_aligned_positions[0]) or (aligned_positions==possible_expected_aligned_positions[1]))

    def test_w_score_extend_gap_without_w(self):
        gap_open=0
        gap_extend=1000
        mrf1 = get_fake_model([0,1,2,3])
        mrf2 = get_fake_model([0,3])
        insert_costs = [
                        {"open":[0]*(mrf1.ncol+1), "extend":[0]*(mrf1.ncol+1)},
                        {"open":[0,gap_open,0], "extend":[0,gap_extend,0]}
                        ]
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models([mrf1,mrf2], self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])
        possible_expected_aligned_positions = [{"pos_ref":[0], "pos_2":[0]}, {"pos_ref":[3], "pos_2":[1]}]
        assert((aligned_positions==possible_expected_aligned_positions[0]) or (aligned_positions==possible_expected_aligned_positions[1]))

    def test_w_score_extend_gap_without_w_mirror(self):
        gap_open=0
        gap_extend=1000
        mrf1 = get_fake_model([0,1,2,3])
        mrf2 = get_fake_model([0,3])
        insert_costs = [
                        {"open":[0,gap_open,0,0,0], "extend":[0,gap_extend,0,0,0]},
                        {"open":[0,gap_open,0], "extend":[0,gap_extend,0]}
                        ]
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models([mrf1,mrf2], self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])


    def test_reverse_w_score_extend_gap_without_w(self):
        gap_open=0
        gap_extend=1000
        mrf1 = get_fake_model([0,3])
        mrf2 = get_fake_model([0,1,2,3])
        insert_costs = [
                        {"open":[0,gap_open,0], "extend":[0,gap_extend,0]},
                        {"open":[0]*(mrf2.ncol+1), "extend":[0]*(mrf2.ncol+1)}
                        ]
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models([mrf1,mrf2], self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(infos_solver['UB'], infos_solver['LB'])
        possible_expected_aligned_positions = [{"pos_ref":[0], "pos_2":[0]}, {"pos_ref":[1], "pos_2":[3]}]
        assert((aligned_positions==possible_expected_aligned_positions[0]) or (aligned_positions==possible_expected_aligned_positions[1]))


    def test_end_gap_free_with_position_specific_insertions(self):
        mrf1 = get_fake_model([0,1])
        mrf2 = get_fake_model([1])
        insert_costs = [
                {"open":[0,0,0], "extend":[0,0,0]},
                {"open":[1000,0], "extend":[0,0]}
                ]
        expected_aligned_positions = {"pos_ref":[1], "pos_2":[0]}
        aligned_positions, aligned_positions_with_gaps, infos_solver = align_two_potts_models([mrf1,mrf2], self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001, free_end_gaps=True)
        self.assertEqual(aligned_positions, expected_aligned_positions)


if __name__=='__main__':
    unittest.main()
