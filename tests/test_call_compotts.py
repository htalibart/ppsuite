import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

from makepotts.potts_object import *
from compotts.call_compotts import *

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


class Test_Call_ComPotts(unittest.TestCase):

    def setUp(self):
        self.output_folder = pathlib.Path(tempfile.mkdtemp())
        self.feature_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.feature_folder)
        self.object = Potts_Object.from_folder(feature_folder=self.feature_folder)
        self.potts_model = self.object.potts_model

    def tearDown(self):
        shutil.rmtree(self.output_folder)
        shutil.rmtree(self.feature_folder)

    def test_align_compotts_object_to_itself(self):
        aligned_positions, infos_solver = align_two_objects([self.object, self.object], self.output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)

    def test_v_score(self):
        aligned_positions, infos_solver = align_two_objects([self.object, self.object], self.output_folder, use_w=False)
        p = self.object.potts_model
        v_score = sum([scalar_product(p.v[i],p.v[i]) for i in range(p.ncol)])
        assert(infos_solver['UB'] - v_score < 1)

    def test_align_small_fake_mrfs(self):
        templates = [["1", "0", "3", "2"],["[0]", "y", "[1]", "[2]"]]
        alnfnames = [self.output_folder/("fake_"+str(i)+".aln") for i in range(2)]
        fastafnames = [self.output_folder/("fake_"+str(i)+".fasta") for i in range(2)]
        crfake.main(templates, alnfnames, fastafnames)
        mrfs = [Potts_Model.from_training_set(fastafnames[i], self.output_folder/("fake_"+str(i)+".mrf")) for i in range(2)]
        aligned_positions, infos_solver = align_two_potts_models(mrfs, self.output_folder)
        self.assertTrue(are_templates_aligned(templates[1], aligned_positions))

    def test_align_different_sizes(self):
        templates = [["y", "y", "y"],["[0]", "y", "[1]", "[2]"]]
        alnfnames = [self.output_folder/("fake_"+str(i)+".aln") for i in range(2)]
        fastafnames = [self.output_folder/("fake_"+str(i)+".fasta") for i in range(2)]
        crfake.main(templates, alnfnames, fastafnames)
        mrfs = [Potts_Model.from_training_set(fastafnames[i], self.output_folder/("fake_diff_size_"+str(i)+".mrf")) for i in range(2)]
        aligned_positions, infos_solver = align_two_potts_models(mrfs, self.output_folder)
        self.assertTrue(are_templates_aligned(templates[1], aligned_positions))
        mrfs.reverse()
        templates.reverse()
        aligned_positions, infos_solver = align_two_potts_models(mrfs, self.output_folder)
        self.assertTrue(are_templates_aligned(templates[1], aligned_positions))

    def test_identity_rescaling_dimensions(self):
        for x in [12, [1, 2, 3], np.ones((3,3))]:
            self.assertTrue(np.array_equal(identity(x), x))

    def test_get_edges_map(self):
        edges_map = get_edges_map(self.potts_model, "none")
        self.assertEqual(edges_map.shape, self.potts_model.w.shape[0:2])
        self.assertEqual(sum([abs(edges_map[i][i]) for i in range(self.potts_model.ncol)]), 0)

    def test_count_edges(self):
        self.assertEqual(count_edges(np.ones((3,3))),9)
        self.assertEqual(count_edges(np.zeros((3,3))),0)


    def test_rescale_mrf(self):
        resc_mrf = get_rescaled_potts_model(self.potts_model, "original_rescaling", "original_rescaling")
        self.assertEqual(self.potts_model.v.shape,resc_mrf.v.shape)
        self.assertEqual(self.potts_model.w.shape,resc_mrf.w.shape)
        self.assertEqual(self.potts_model.ncol, resc_mrf.ncol)


    def test_align_rescaled_mrf_to_itself(self):
        resc_mrf = get_rescaled_potts_model(self.potts_model, "original_rescaling", "original_rescaling")
        aligned_positions, infos_solver = align_two_potts_models([resc_mrf, resc_mrf], self.output_folder)
        similarity_global = infos_solver["similarity_global"]


    def test_alpha_w(self):
        alpha_w = 2
        aligned_positions_alpha, infos_solver_alpha = align_two_potts_models([self.potts_model, self.potts_model], self.output_folder, alpha_w=alpha_w, use_v=False)
        aligned_positions_normal, infos_solver_normal = align_two_potts_models([self.potts_model, self.potts_model], self.output_folder, use_v=False)
        assert((infos_solver_alpha['UB']-alpha_w*infos_solver_normal['UB'])<1)



if __name__=='__main__':
    unittest.main()
