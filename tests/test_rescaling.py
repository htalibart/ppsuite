import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils import files_management as fm
from makepotts.potts_object import *
from makepotts.rescaling import *
from makepotts.potts_model import *

class Test_Rescaling(unittest.TestCase):

    def setUp(self):
        self.feature_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.feature_folder)
        self.potts_model = Potts_Model.from_msgpack(self.feature_folder/"potts_model.mrf")


    def tearDown(self):
        #shutil.rmtree(self.feature_folder)
        pass

    def test_rescaling(self):
        cf = Potts_Object.from_files(self.feature_folder, aln_file=SMALL_ALN_1CC8, trim_alignment=False, max_nb_sequences=250, v_rescaling_function="original_rescaling", w_rescaling_function="original_rescaling") 
       
    def test_simulate_uniform_pc_on_v(self):
        tau = 0.2
        resc = get_rescaled_potts_model(self.potts_model, "simulate_uniform_pc_on_v", "identity", use_w=True, v_rescaling_tau=tau)
        assert(self.potts_model.v.shape==resc.v.shape)


    def test_simulate_uniform_pc_on_w(self):
        tau = 0.2
        resc = get_rescaled_potts_model(self.potts_model, "identity", "simulate_uniform_pc_on_w", use_w=True, w_rescaling_tau=0.5, beta_softmax_w=10, w_back_to_scale=True)


    def test_simulate_uniform_pc_on_w_many_parameters(self):
        tau = 0.2
        str_="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwyz"
        my_kwargs = {str_[i]:i for i in range(len(str_))}
        resc = get_rescaled_potts_model(self.potts_model, "identity", "simulate_uniform_pc_on_w", use_w=True, w_rescaling_tau=0.5, beta_softmax_w=10, w_back_to_scale=True, **my_kwargs)

if __name__=='__main__':
    unittest.main()
