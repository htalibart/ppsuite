from comutils.util import *
import numpy as np
import math
from numpy import linalg as LA

# TODO vectoriser ou edges_map
# TODO décider d'un seuil
# TODO commenter

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


def compute_v_scores(mrf1, mrf2, v_score_function, v_coeff=1, offset_v=0, **kwargs):
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = v_score_function(mrf1.v[i], mrf2.v[k])-offset_v
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


def get_w_score_at(w_scores,i,k,j,l):
    if (i==j):
        return 0
    elif (k==l):
        return 0
    else:
        if (l<k):
            if (j<i):
                return w_scores[j+int(i*(i+1)/2)][l+int(k*(k+1)/2)]
            else:
                return get_w_score_at(w_scores,j,k,i,l)
        else:
            return get_w_score_at(w_scores,i,l,j,k)



def get_epsilon(precision_method, selfscores):
    if precision_method.startswith("arbitrary_"):
        return float(precision_method[len("arbitrary_"):])
    elif precision_method.startswith("similarity_"):
        sim_diff = float(precision_method[len("similarity_"):])
        return (sim_diff/2)*sum(selfscores)
    else:
        raise Exception("Undefined epsilon strategy")



#def compute_scores_etc(mrfs, v_score_function=scalar_product, w_score_function=scalar_product, use_v=True, use_w=True, alpha_w=1, w_threshold_method="percentile_0", precision_method="similarity_0.005", **kwargs):
#    if use_w:
#        edges_maps = [get_edges_map(mrf, w_threshold_method) for mrf in mrfs]
#    else:
#        edges_maps = [np.zeros((mrf.w.shape[0:2])) for mrf in mrfs]
#    if use_v:
#        v_scores = compute_v_scores(*mrfs, v_score_function, v_coeff=1)
#    else:
#        v_scores = np.zeros(tuple([mrf.v.shape[0] for mrf in mrfs]))
#    if use_w:
#        w_scores = compute_w_scores(*mrfs, *edges_maps, w_score_function, w_coeff=alpha_w)
#    else:
#        w_scores = np.zeros(1)
#    selfscores = [compute_selfscore(mrf, edges_map, v_score_function, w_score_function, use_v, use_w, alpha_w, **kwargs) for mrf, edges_map in zip(mrfs, edges_maps)]
#    epsilon = get_epsilon(precision_method, selfscores)
#    return v_scores, w_scores, edges_maps, selfscores, epsilon
#

def compute_self_w_scores(mrf, edges_map, w_score_function, **kwargs):
    w_score = 0
    for i in range(mrf.ncol-1):
        for j in range(i+1,mrf.ncol):
            if edges_map[i][j]:
                w_score+=w_score_function(mrf.w[i][j],mrf.w[i][j])
    return w_score


def compute_selfscore(mrf, edges_map, use_v=True, use_w=True, alpha_w=1, offset_v=0, **kwargs):
    v_score_function = scalar_product
    w_score_function = scalar_product
    if use_v:
        v_score = sum([v_score_function(vi,vi)-offset_v for vi in mrf.v])
    else:
        v_score = 0
    if use_w:
        w_score = compute_self_w_scores(mrf, edges_map, w_score_function)
    else:
        w_score = 0
    selfcomp = v_score+alpha_w*w_score
    return selfcomp
