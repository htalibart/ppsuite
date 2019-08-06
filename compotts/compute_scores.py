from basic_modules.util import *
import numpy as np
import math
from numpy import linalg as LA

# TODO vectoriser ou edges_map
# TODO décider d'un seuil

def count_edges(edges_map):
    return sum(edges_map.flatten())


def get_w_threshold(mrf, w_threshold_method):
    if w_threshold_method.startswith("percentile_"):
        p = int(w_threshold_method[len("percentile_"):])
        return np.percentile(mrf.get_w_norms(), p)
    elif w_threshold_method=="magnitude":
        w_max = np.max(mrf.get_w_norms().flatten())
        return w_max/10
    else:
        return 0


def get_edges_map(mrf, w_threshold_method):
    w_threshold = get_w_threshold(mrf, w_threshold_method)
    return 1*(mrf.get_w_norms()>w_threshold)


def compute_v_scores(mrf1, mrf2, v_score_function, v_coeff=1, **kwargs):
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = v_score_function(mrf1.v[i], mrf2.v[k])
    return v_coeff*v_scores

def compute_w_scores(mrf1, mrf2, edges_map1, edges_map2, w_score_function, w_coeff=1, **kwargs):
    """ symmetric matrix : w[i][j]=v[j+i*(i+1)/2])"""
    len1 = int(mrf1.ncol*(mrf1.ncol+1)/2)
    len2 = int(mrf2.ncol*(mrf2.ncol+1)/2)
    print("computing w scores (LA="+str(mrf1.ncol)+", LB="+str(mrf2.ncol)+" : "+str(len1)+"x"+str(len2)+")")
    w_scores = np.zeros((len1,len2), np.dtype('f4'))
    for i in range(mrf1.ncol):
        for j in range(i+1):
            if edges_map1[i][j]:
                for k in range(mrf2.ncol):
                    for l in range(k+1):
                        if edges_map2[k][l]:
                            w_scores[j+int(i*(i+1)/2)][l+int(k*(k+1)/2)] = w_score_function(mrf1.w[i][j], mrf2.w[k][l])
    return w_coeff*w_scores


def get_vw_coeffs(mrfs, vw_coeff_method, edges_maps, v_score_function=scalar_product, w_score_function=scalar_product, use_v=True, use_w=True):
    if vw_coeff_method.startswith("arbitrary_"):
        return [float(strcoeff) for strcoeff in vw_coeff_method[len("arbitrary_"):].split('_')]
    elif vw_coeff_method.startswith("scoremax_"):
        hmax = float(vw_coeff_method[len("scoremax_"):])
        w_norms = [mrf.get_w_norms()*edge_map for mrf, edge_map in zip(mrfs, edges_maps)]
        denom = mrfs[0].get_v_norm()*mrfs[1].get_v_norm()+(1/2)*LA.norm(w_norms[0])*LA.norm(w_norms[1])
        alpha = hmax/denom
        #alpha = 2*hmax/sum([compute_selfscore(mrf, edges_map, v_score_function, w_score_function, use_v, use_w, v_coeff=1, w_coeff=1) for mrf, edges_map in zip(mrfs, edges_maps)])
        return [alpha, alpha]
    else:
        return [1, 1]



def get_gap_costs(gap_cost_method, v_scores, vw_coeff_method):
    if gap_cost_method.startswith("arbitrary_"):
        return [float(strcoeff) for strcoeff in gap_cost_method[len("arbitrary_"):].split('_')]
    elif gap_cost_method.startswith("one_scoremax_pos") and (vw_coeff_method.startswith("scoremax_")):
        hmax = float(vw_coeff_method[len("scoremax_"):])
        gap_open = hmax/min(v_scores.shape[0], v_scores.shape[1])
        return [gap_open,0]
    elif gap_cost_method.startswith("max_score_v_times_"):
        fact = float(gap_cost_method[len("max_score_v_times_"):])
        gap_open = fact*np.max(v_scores.flatten())
        return[float(gap_open), 0]
    else:
        print("No gap cost method, returning [5,0]")
        return [5,0]


def get_epsilon(epsilon_method, vw_coeff_method):
    if epsilon_method.startswith("arbitrary_"):
        return float(epsilon_method[len("arbitrary_"):])
    elif (vw_coeff_method.startswith("scoremax_")) and (epsilon_method.startswith("p_scoremax_")):
        Hmax = float(vw_coeff_method[len("scoremax_"):])
        p = float(epsilon_method[len("p_scoremax_"):])
        return p*Hmax
    else:
        print("No precision method, returning 1")
        return 1



def compute_scores_etc(mrfs, v_score_function=scalar_product, w_score_function=scalar_product, use_v=True, use_w=True, vw_coeff_method="arbitrary_1_1", w_threshold_method="percentile_0",gap_cost_method="arbitrary_5_0", epsilon_method="arbitrary_1", **kwargs):
    if use_w:
        edges_maps = [get_edges_map(mrf, w_threshold_method) for mrf in mrfs]
    else:
        edges_maps = [np.zeros((mrf.w.shape[0:2])) for mrf in mrfs]
    [v_coeff, w_coeff] = get_vw_coeffs(mrfs, vw_coeff_method, edges_maps, v_score_function, w_score_function, use_v, use_w)
    if use_v:
        v_scores = compute_v_scores(*mrfs, v_score_function, v_coeff=v_coeff)
    else:
        v_scores = np.zeros(tuple([mrf.v.shape[0] for mrf in mrfs]))
    if use_w:
        w_scores = compute_w_scores(*mrfs, *edges_maps, w_score_function, w_coeff=w_coeff)
    else:
        w_scores = np.zeros(1)
    [gap_open, gap_extend] = get_gap_costs(gap_cost_method, v_scores, vw_coeff_method)
    selfscores = [compute_selfscore(mrf, edges_map, v_score_function, w_score_function, use_v, use_w, vw_coeff_method, **kwargs) for mrf, edges_map in zip(mrfs, edges_maps)]
    epsilon = get_epsilon(epsilon_method, vw_coeff_method)
    return v_scores, w_scores, edges_maps, selfscores, gap_open, gap_extend, epsilon


def compute_self_w_scores(mrf, edges_map, w_score_function, **kwargs):
    w_score = 0
    for i in range(mrf.ncol-1):
        for j in range(i+1,mrf.ncol):
            if edges_map[i][j]:
                w_score+=w_score_function(mrf.w[i][j],mrf.w[i][j])
    return w_score


def compute_selfscore(mrf, edges_map, v_score_function, w_score_function, use_v=True, use_w=True, vw_coeff_method="arbitrary_1_1", **kwargs):
    [v_coeff, w_coeff] = get_vw_coeffs([mrf,mrf], vw_coeff_method, [edges_map, edges_map], v_score_function, w_score_function, use_v, use_w)
    if use_v:
        v_score = sum([v_score_function(vi,vi) for vi in mrf.v])
    else:
        v_score = 0
    if use_w:
        w_score = compute_self_w_scores(mrf, edges_map, w_score_function)
    else:
        w_score = 0
    selfcomp = v_coeff*v_score+w_coeff*w_score
    return selfcomp
