import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils import files_management as fm
from makepotts.potts_object import *
from makepotts.rescaling import *
from makepotts.potts_model import *

import time

class Test_Rescaling_CPP(unittest.TestCase):

    def setUp(self):
        self.potts_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.potts_folder)
        self.potts_model = Potts_Model.from_msgpack(self.potts_folder/"potts_model.mrf")


    def tearDown(self):
        shutil.rmtree(self.potts_folder)
        pass

    def test_rescale_w_cpp(self):
        w = self.potts_model.w

        initial_w = w.copy()
        start = time.time()
        rescaled_w_cpp = simulate_uniform_pc_on_w_with_cpp(w)
        end = time.time()
        time_cpp = end-start
        assert(np.array_equal(w, initial_w))

        start = time.time()
        rescaled_w_python = simulate_uniform_pc_on_w_with_python(w)
        end = time.time()
        time_python = end-start

        print("time C++:", time_cpp, "time Python:", time_python)

        assert(np.allclose(rescaled_w_cpp, rescaled_w_python, atol=1e-5))


    def test_rescale_v_cpp(self):
        v = self.potts_model.v

        initial_v = v.copy()
        start = time.time()
        rescaled_v_cpp = simulate_uniform_pc_on_v_with_cpp(v)
        end = time.time()
        time_cpp = end-start
        assert(np.array_equal(v, initial_v))

        start = time.time()
        rescaled_v_python = simulate_uniform_pc_on_v_with_python(v)
        end = time.time()
        time_python = end-start

        print("time C++:", time_cpp, "time Python:", time_python)
    
        assert(np.allclose(rescaled_v_cpp, rescaled_v_python, atol=1e-5))


if __name__=='__main__':
    unittest.main()
