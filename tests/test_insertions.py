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

def get_vi_conserved_letter(conserved_letter):
    via_conserved = 1
    vi = np.ones(21)*(-via_conserved/(21-1))
    vi[conserved_letter] = via_conserved
    return vi


def get_fake_model(conserved_letters_list):
    v = np.array([get_vi_conserved_letter(c) for c in conserved_letters_list])
    w = np.zeros((len(conserved_letters_list), len(conserved_letters_list), 21, 21))
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
        mrfs = [get_fake_model([0,1,2]), get_fake_model([3,4,1,2])]
        insert_costs = [{"open":[0,1000,0,0],"extend":[0,0,0,0]}, {"open":[0,0,0,0,0],"extend":[0,0,0,0,0]}]
        expected_aligned_positions = {"pos_ref":[1,2], "pos_2":[2,3]}
        aligned_positions, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(aligned_positions, expected_aligned_positions)


    def test_costly_extend(self):
        mrfs = [get_fake_model([0,1,2,3]), get_fake_model([0,4,5,6,1,2,3])]
        insert_costs = [{"open":[0,1,0,0,0],"extend":[0,0.5,0,0,0]}, {"open":[0,0,0,0,0,0,0,0],"extend":[0,0,0,0,0,0,0,0]}]
        expected_aligned_positions = {"pos_ref":[1,2,3], "pos_2":[4,5,6]}
        aligned_positions, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
        self.assertEqual(aligned_positions, expected_aligned_positions)





#    def test_costly_open_with_templates(self):
#        templates = [['y','y','y'],['[0]', 'y', '[1]', '[2]']]
#        insert_costs = [{"open":[0,1000,0,0],"extend":[0,0,0,0]}, {"open":[0,0,0,0,0],"extend":[0,0,0,0,0]}]
#        expected_aligned_positions = {"pos_ref":[1,2], "pos_2":[2,3]}
#        alnfnames = [self.output_folder/("fake_"+str(i)+".aln") for i in range(2)]
#        fastafnames = [self.output_folder/("fake_"+str(i)+".fasta") for i in range(2)]
#        crfake.main(templates, alnfnames, fastafnames, n=500)
#        mrfs = [Potts_Model.from_training_set(fastafnames[i], self.output_folder/("fake_"+str(i)+".mrf")) for i in range(2)]
#        aligned_positions, infos_solver = align_two_potts_models(mrfs, self.output_folder, insert_costs=insert_costs, sim_min=-100, epsilon_sim=0.0001)
#        self.assertEqual(aligned_positions, expected_aligned_positions)
#

if __name__=='__main__':
    unittest.main()
