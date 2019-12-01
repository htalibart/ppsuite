import unittest
import shutil, tempfile
import pathlib

import numpy as np

from comfeature.comfeature import *
from compotts.manage_positions import *

import pkg_resources
from tests.resources_manager import *

class Test_ManagePositions(unittest.TestCase):

    def setUp(self):
        self.output_folder = pathlib.Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.output_folder)

    def test_get_alignment_with_gaps(self):
        aligned_positions = {"pos_ref":[4,5,7,8,9,10], "pos_2":[0,1,2,3,4,5]}
        alignment_with_gaps = {"pos_ref":[0,1,2,3,4,5,6,7,8,9,10], "pos_2":['-','-','-','-',0,1,'-',2,3,4,5]}
        res = get_alignment_with_gaps(aligned_positions)
        self.assertEqual(alignment_with_gaps, res)

    def test_get_alignment_with_gaps_2(self):
        aligned_positions = {"pos_ref":[4,5,8,9,10], "pos_2":[0,1,3,4,5]}
        alignment_with_gaps = {"pos_ref":[0,1,2,3,4,5,'X',8,9,10], "pos_2":['-','-','-','-',0,1,'X',3,4,5]}
        res = get_alignment_with_gaps(aligned_positions)
        self.assertEqual(alignment_with_gaps, res)

    def test_get_alignment_with_gaps_2(self):
        aligned_positions = {"pos_ref":[0,1,2,3], "pos_2":[0,2,3,4]}
        alignment_with_gaps = {"pos_ref":[0,'-',1,2,3], "pos_2":[0,1,2,3,4]}
        res = get_alignment_with_gaps(aligned_positions)
        self.assertEqual(alignment_with_gaps, res)

    def test_get_alignment_with_gaps_3(self):
        aligned_positions = {"pos_ref":[4,5,8,9,10], "pos_2":[1,2,3,4,5]}
        alignment_with_gaps = {"pos_ref":['X',4,5,6,7,8,9,10], "pos_2":['X',1,2,'-','-',3,4,5]}
        res = get_alignment_with_gaps(aligned_positions)
        self.assertEqual(alignment_with_gaps, res)

    def test_get_seqs_aligned(self):
        objs = []
        for k in range(2):
            seq_file = FAKE_SEQS_FOLDER/("fake_seq_"+str(k)+".fasta")
            feature_folder = pathlib.Path(tempfile.mkdtemp())
            objs.append(ComFeature.from_files(feature_folder=feature_folder, sequence_file=seq_file, inference_type="one_hot"))
        aligned_positions = {"pos_ref":[4,5,7,8,9,10], "pos_2":[0,1,2,3,4,5]}
        seqs_aligned = get_seqs_aligned(aligned_positions, objs)
        self.assertEqual(seqs_aligned[0], 'YFYFMAEIKEH')
        self.assertEqual(seqs_aligned[1], '----MA-IKDH')



if __name__=='__main__':
    unittest.main()
