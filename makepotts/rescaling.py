from math import *
from comutils.global_variables import *
import numpy as np
from comutils.util import *
from makepotts.potts_model import *

VECTORIZABLE_FUNCTIONS = ["identity", "original_rescaling", "symmetric_relu_like", "shifted_relu", "threshold_on_wijab", "exponential"]


def get_rescaled_potts_model(potts_model, v_rescaling_function_name, w_rescaling_function_name, use_w=True, **kwargs):
    """ returns a copy of Potts Model @mrf rescaled using @v_rescaled_function_name rescaling function to rescale v and @w_rescaled_function for w if @use_w """
    useful_kwargs = {}
    for key in kwargs.keys():
        if key in ["alpha_rescaling", "wijab_threshold", "v_rescaling_tau", "w_rescaling_tau", "v_back_to_scale", "w_back_to_scale", "beta_softmax_w"]:
            useful_kwargs[key] = kwargs[key]
    if (v_rescaling_function_name=="identity") and (w_rescaling_function_name=="identity"):
        return potts_model
    else:
        print("rescaling Potts model")
        t_v = get_rescaled_parameters(potts_model.v, v_rescaling_function_name, **useful_kwargs)
        t_v[:,20]=0 # keep gap parameter to 0
        if use_w:
            t_w = get_rescaled_parameters(potts_model.w, w_rescaling_function_name, **useful_kwargs)
        else:
            t_w = np.zeros_like(potts_model.w)
        return Potts_Model.from_parameters(t_v, t_w, name=potts_model.name+'_'+v_rescaling_function_name+'_'+w_rescaling_function_name)


def get_rescaled_parameters(x, rescaling_function_name, **kwargs):
    if rescaling_function_name in VECTORIZABLE_FUNCTIONS:
        vfunc = np.vectorize(eval(rescaling_function_name))
        return vfunc(x, **kwargs)
    else:
        return eval(rescaling_function_name)(x,**kwargs)


def get_background_v0():
    f0 = np.array(AA_BACKGROUND_FREQUENCIES)
    logf0 = np.log(f0)
    vaa = logf0-(1/len(f0))*np.sum(logf0)
    v0 = np.append(vaa, [0])
    return v0

def remove_v0(x, v_rescaling_function_name="identity", **kwargs):
    v0 = get_background_v0()
    v0_rescaled = get_rescaled_parameters(v0, v_rescaling_function_name, **kwargs)
    return x-v0_rescaled


def identity(x, **kwargs):
    return x

def original_rescaling(x, alpha_rescaling=1, **kwargs):
    return sign_ind(x)*(exp(alpha_rescaling*abs(x))-1)


def threshold_on_wijab(x, wijab_threshold=0.05, **kwargs):
        return x*(abs(x)>=wijab_threshold)

def exponential(x, **kwargs):
    return np.exp(x)


def simulate_uniform_pc_on_v(v, v_rescaling_tau=1/2, v_back_to_scale=False, **kwargs):
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
            resc_tmp[a] = log((1-v_rescaling_tau)*exp(vi[a])/S + v_rescaling_tau/q)

        resc_vi = np.zeros_like(vi)
        S_all = np.sum(resc_tmp)
        for a in range(q):
            resc_vi[a] = resc_tmp[a]-(1/q)*S_all
        if v_back_to_scale and euclidean_norm(resc_vi)!=0:
            resc_vi = euclidean_norm(vi)/euclidean_norm(resc_vi)*resc_vi
        resc_v[i] = resc_vi
    return resc_v


def simulate_uniform_pc_on_wij(w, rescaling_tau=0.5, beta=10, w_back_to_scale=False, **kwargs):
    w_flat = w.flatten()
    S = sum([exp(beta*elt) for elt in w_flat])
    
    resc_tmp = np.zeros_like(w_flat)
    for elt in range(len(w_flat)):
        frac = exp(beta*w_flat[elt])/S
        resc_tmp[elt] = log((1-rescaling_tau)*frac+rescaling_tau/len(w_flat))
    
    resc_flat = np.zeros_like(w_flat)
    S_all = np.sum(resc_tmp)
    for elt in range(len(w_flat)):
        resc_flat[elt] = (1/beta)*(resc_tmp[elt]-(1/len(w_flat))*S_all)

    resc_unflat = resc_flat.reshape(w.shape)
    
    if w_back_to_scale and euclidean_norm(resc_unflat)!=0:
        resc_unflat = euclidean_norm(w)/euclidean_norm(resc_unflat)*resc_unflat
    return resc_unflat


def simulate_uniform_pc_on_w(w, w_rescaling_tau=0.5, beta_softmax_w=10, w_back_to_scale=False, **kwargs):
    resc_w = np.zeros_like(w)
    for i in range(len(resc_w)):
        for j in range(len(resc_w)):
            resc_w[i,j] = simulate_uniform_pc_on_wij(w[i][j], rescaling_tau=w_rescaling_tau, beta=beta_softmax_w,
                                                            w_back_to_scale=w_back_to_scale)
    return resc_w
