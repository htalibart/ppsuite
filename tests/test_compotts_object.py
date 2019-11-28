import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from compotts.compotts_object import *


class Test_ComPotts_Object(unittest.TestCase):

    def test_from_mrf(self):
        co = ComPotts_Object(potts_model_file=MRF_1CC8)
        shutil.rmtree(co.folder)

    def test_from_folder_with_everything(self):
        input_folder_name = '/tmp/'+next(tempfile._get_candidate_names())
        input_folder = pathlib.Path(input_folder_name)
        shutil.copytree(RESOURCES_1CC8_EVERYTHING_FOLDER, input_folder)
        co = ComPotts_Object(input_folder=input_folder)
        shutil.rmtree(input_folder)

    def test_from_nothing_with_seq_and_a3m(self):
        input_folder = pathlib.Path(tempfile.mkdtemp())
        co = ComPotts_Object(seq_file=SEQ_1CC8, a3m_file=A3M_1CC8, input_folder=input_folder)
        shutil.rmtree(co.folder)

    def test_from_folder_with_seq_and_a3m(self):
        input_folder = pathlib.Path(tempfile.mkdtemp())
        shutil.copy(SEQ_1CC8, input_folder)
        shutil.copy(A3M_1CC8, input_folder)
        co = ComPotts_Object(input_folder=input_folder)
        shutil.rmtree(input_folder)



if __name__=='__main__':
    unittest.main()
