from math import *
import numpy as np
from comutils.util import *
from comutils.potts_model import *

def get_rescaled_potts_model(mrf, rescaling_function_name, use_w=True, **kwargs):
    """ returns a copy of Potts Model @mrf rescaled using @rescaled_function_name rescaling function, and using or not w """
    rescaling_function = eval(rescaling_function_name)
    t_v = np.zeros_like(mrf.v)
    for i in range(len(mrf.v)):
        t_v[i] = rescale_parameter(mrf.v[i], rescaling_function, parameter_type="v", **kwargs)
    t_w = np.zeros_like(mrf.w)
    if use_w:
        for i in range(len(mrf.w)):
            for j in range(len(mrf.w)):
                t_w[i][j] = rescale_parameter(mrf.w[i][j], rescaling_function, parameter_type="w", **kwargs)
    return Potts_Model.from_parameters(t_v, t_w, name=mrf.name+"_"+rescaling_function_name)


def rescale_parameter(x, rescaling_function, **kwargs):
    """ rescale one parameter with @rescaling_function """
    vfunc = np.vectorize(rescaling_function)
    return vfunc(x, **kwargs)



# AVAILABLE RESCALING FUNCTIONS

def identity(x, **kwargs):
    return x

def original_rescaling(x, **kwargs):
    return sign_ind(x)*(exp(abs(x))-1)


def symmetric_relu_like(x, v_lower_threshold=-2, v_upper_threshold=3, w_lower_threshold=-0.01, w_upper_threshold=0.01, **kwargs):
    if kwargs["parameter_type"]=="v":
        lower_threshold = v_lower_threshold
        upper_threshold = v_upper_threshold
    else :
        lower_threshold = w_lower_threshold
        upper_threshold = w_upper_threshold
    return x*((x>=upper_threshold) or (x<=lower_threshold))


def shifted_relu(x, v_threshold=1, w_threshold=0.01, **kwargs):
    if kwargs["parameter_type"]=="v":
        threshold = v_threshold
    else:
        threshold = w_threshold
    return x*(x>=threshold)

def add_number(x, shift=3, **kwargs):
    if kwargs["parameter_type"]=="v":
        return x+shift
    else:
        return x
