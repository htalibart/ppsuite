from math import *
import numpy as np
from comutils.util import *
from makepotts.potts_model import *

VECTORIZABLE_FUNCTIONS = ["identity", "original_rescaling", "symmetric_relu_like", "shifted_relu", "add_number", "threshold_on_wijab"]


def get_rescaled_potts_model(potts_model, v_rescaling_function_name, w_rescaling_function_name, use_w=True, **kwargs):
    """ returns a copy of Potts Model @mrf rescaled using @v_rescaled_function_name rescaling function to rescale v and @w_rescaled_function for w if @use_w """
    if (v_rescaling_function_name=="identity") and (w_rescaling_function_name=="identity"):
        return potts_model
    else:
        print("rescaling Potts model")
        t_v = get_rescaled_parameters(potts_model.v, v_rescaling_function_name, **kwargs)
        t_v[:,20]=0 # keep gap parameter to 0
        if use_w:
            t_w = get_rescaled_parameters(potts_model.w, w_rescaling_function_name, **kwargs)
        else:
            t_w = np.zeros_like(potts_model.w)
        return Potts_Model.from_parameters(t_v, t_w, name=potts_model.name+'_'+v_rescaling_function_name+'_'+w_rescaling_function_name)


def get_rescaled_parameters(x, rescaling_function_name, **kwargs):
    if rescaling_function_name in VECTORIZABLE_FUNCTIONS:
        vfunc = np.vectorize(eval(rescaling_function_name))
        return vfunc(x, **kwargs)
    else:
        return eval(rescaling_function_name)(x,**kwargs)


def identity(x, **kwargs):
    return x

def original_rescaling(x, alpha_rescaling=1, **kwargs):
    return sign_ind(x)*(exp(alpha_rescaling*abs(x))-1)

def add_number(x, v_shift=3, **kwargs):
        return x+v_shift

def threshold_on_wijab(x, wijab_threshold=0.05, **kwargs):
        return x*(abs(x)>=wijab_threshold)

def simulate_uniform_pc_on_v(v, rescaling_tau=1/2, **kwargs):
    resc_v = np.zeros_like(v)
    for i in range(len(v)):
        vi = v[i]
        resc_vi = np.zeros_like(vi)
        q=20
        S = 0
        for b in range(q):
            S+=exp(vi[b])

        resc_tmp= np.zeros_like(vi)
        for a in range(q):
            resc_tmp[a] = log((1-rescaling_tau)*exp(vi[a])/S + rescaling_tau/q)

        resc_vi = np.zeros_like(vi)
        S_all = np.sum(resc_tmp)
        for a in range(q):
            resc_vi[a] = resc_tmp[a]-(1/q)*S_all
        resc_v[i] = resc_vi
    return resc_v
