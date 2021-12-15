import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from makepotts.potts_object import *
from comutils import files_management as fm

class Test_MakePotts_Main(unittest.TestCase):

    def setUp(self):
        self.potts_folder = pathlib.Path(tempfile.mkdtemp())
        pass

    def tearDown(self):
        shutil.rmtree(self.potts_folder)
        pass


    def test_from_hhblits_files_a3m_only(self):
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--aln_with_insertions", str(A3M_1CC8)]
        potts_object = main(makepotts_args)
        assert(potts_object.potts_model_file.is_file())


    def test_from_aln_file(self):
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--aln_file", str(ALN_1CC8)]
        potts_object = main(makepotts_args)
        assert(potts_object.potts_model_file.is_file())


    def test_use_less_sequences_and_dont_infer(self):
        nb_sequences=10
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--aln_with_insertions", str(A3M_1CC8), "--max_nb_sequences", str(nb_sequences), "--dont_infer_potts_model"]
        potts_object = main(makepotts_args)
        assert(not potts_object.potts_model_file.is_file())
        assert(fm.get_nb_sequences_in_fasta_file(potts_object.aln_train)==nb_sequences)


    def test_one_hot(self):
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--inference_type", "one_hot"]
        potts_object = main(makepotts_args)
        assert(potts_object.potts_model_file.is_file())

    def test_one_submat(self):
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--inference_type", "one_submat"]
        potts_object = main(makepotts_args)
        assert(potts_object.potts_model_file.is_file())

    def test_light(self):
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--aln_with_insertions", str(A3M_1CC8), "--light"]
        potts_object = main(makepotts_args)
        assert(potts_object.potts_model_file.is_file())


    def test_infer_mfdca(self):
        makepotts_args = ["--potts_folder", str(self.potts_folder), "--sequence_file", str(SEQ_1CC8), "--aln_with_insertions", str(A3M_1CC8), "--inference_method", "mfDCA"]
        potts_object = main(makepotts_args)
        assert(potts_object.potts_model_file.is_file())

if __name__=='__main__':
    unittest.main()
