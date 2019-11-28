import unittest
import shutil, tempfile
import os
import pathlib
import numpy as np

from compotts.compotts_object import *
from compotts.call_compotts import *

import pkg_resources
from tests.resources_manager import *

import basic_modules.create_fake_data as crfake


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
        self.input_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(RESOURCES_1CC8_EVERYTHING_FOLDER, self.input_folder)
        self.object = ComPotts_Object(input_folder=self.input_folder)
        self.mrf = self.object.mrf

    def tearDown(self):
        shutil.rmtree(self.output_folder)
        shutil.rmtree(self.input_folder)

    def test_align_compotts_object_to_itself(self):
        aligned_positions, infos_solver = align_two_objects([self.object, self.object], self.output_folder)
        similarity_global = infos_solver["similarity_global"]
        self.assertEqual(similarity_global,1)

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
        edges_map = get_edges_map(self.mrf, "none")
        self.assertEqual(edges_map.shape, self.mrf.w.shape[0:2])
        self.assertEqual(sum([abs(edges_map[i][i]) for i in range(self.mrf.ncol)]), 0)

    def test_count_edges(self):
        self.assertEqual(count_edges(np.ones((3,3))),9)
        self.assertEqual(count_edges(np.zeros((3,3))),0)


    def test_rescale_mrf(self):
        resc_mrf = get_rescaled_mrf(self.mrf, "original_rescaling")
        self.assertEqual(self.mrf.v.shape,resc_mrf.v.shape)
        self.assertEqual(self.mrf.w.shape,resc_mrf.w.shape)
        self.assertEqual(self.mrf.ncol, resc_mrf.ncol)


    def test_align_rescaled_mrf_to_itself(self):
        resc_mrf = get_rescaled_mrf(self.mrf, "original_rescaling")
        aligned_positions, infos_solver = align_two_potts_models([resc_mrf, resc_mrf], self.output_folder)
        similarity_global = infos_solver["similarity_global"]

if __name__=='__main__':
    unittest.main()
