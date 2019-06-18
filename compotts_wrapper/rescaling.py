from math import *
import numpy as np
from util import *
from potts_model import *

def original_rescaling(x, **kwargs):# TODO test
    return sign_ind(x)*(exp(abs(x))-1)

def get_rescaled_mrf(mrf, rescaling_function_name, **kwargs): # TODO test
    rescaling_function = eval(rescaling_function_name)
    t_v = np.zeros_like(mrf.v)
    for i in range(len(mrf.v)):
        t_v[i] = rescale_parameter(mrf.v[i], rescaling_function)
    t_w = np.zeros_like(mrf.w)
    for i in range(len(mrf.w)):
        for j in range(len(mrf.w)):
            t_w[i][j] = rescale_parameter(mrf.w[i][j], rescaling_function)
    return Potts_Model.from_parameters(t_v, t_w, name=mrf.name+"_"+rescaling_function_name)

def rescale_parameter(x, rescaling_function, **kwargs): # TODO test
    vfunc = np.vectorize(rescaling_function)
    return vfunc(x, **kwargs)


def identity(x, **kwargs):
    return x

