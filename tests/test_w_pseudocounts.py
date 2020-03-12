import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils import files_management as fm
from makepotts.potts_object import *
from makepotts.potts_model import *
from makepotts.pseudo_w import *

class Test_Pseudo_W(unittest.TestCase):

    def setUp(self):
        self.feature_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.feature_folder)
        self.potts_model = Potts_Model.from_msgpack(self.feature_folder/"potts_model.mrf")


    def tearDown(self):
        #shutil.rmtree(self.feature_folder)
        pass

       
    def test_simulate_submat_on_w(self):
        tau = 0.5
        resc = get_potts_model_with_pseudo_w(self.potts_model, tau)


    
    def test_potts_model_with_pseudo_w_and_rescale(self):
        w_submat_tau=0.001
        resc = get_potts_model_with_pseudo_w_and_rescale(self.potts_model, w_submat_tau)
        resc.to_msgpack("tmp_test.mrf")

if __name__=='__main__':
    unittest.main()
