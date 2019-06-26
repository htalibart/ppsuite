from basic_modules.util import *
import numpy as np
import math

# TODO vectoriser ou edges_map
# TODO décider d'un seuil

def count_edges(edges_map):
    return sum(edges_map.flatten())


def get_w_threshold(mrf, w_threshold_method):
    if w_threshold_method.startswith("percentile_"):
        p = int(w_threshold_method[len("percentile_"):])
        return np.percentile(mrf.get_w_norms(), p)
    else:
        return 0


def get_edges_map(mrf, w_threshold_method):
    w_threshold = get_w_threshold(mrf, w_threshold_method)
    return 1*(mrf.get_w_norms()>w_threshold)


def compute_v_scores(mrf1, mrf2, v_score_function, **kwargs):
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = v_score_function(mrf1.v[i], mrf2.v[k])
    return v_scores


def compute_w_scores(mrf1, mrf2, edges_map1, edges_map2, w_score_function, **kwargs):
    """ symmetric matrix : w[i][j]=v[j+i*(i+1)/2])"""
    len1 = int(mrf1.ncol*(mrf1.ncol+1)/2)
    len2 = int(mrf2.ncol*(mrf2.ncol+1)/2)
    print("computing w scores (LA="+str(mrf1.ncol)+", LB="+str(mrf2.ncol)+" : "+str(len1)+"x"+str(len2)+")")
    w_scores = np.zeros((len1,len2))
    for i in range(mrf1.ncol):
        for j in range(i+1):
            if edges_map1[i][j]:
                for k in range(mrf2.ncol):
                    for l in range(k+1):
                        if edges_map2[k][l]:
                            w_scores[j+int(i*(i+1)/2)][l+int(k*(k+1)/2)] = w_score_function(mrf1.w[i][j], mrf2.w[k][l])
    return w_scores


def compute_self_w_scores(mrf, edges_map, w_score_function, **kwargs):
    w_score = 0
    for i in range(mrf.ncol-1):
        for j in range(i+1,mrf.ncol):
            if edges_map[i][j]:
                w_score+=w_score_function(mrf.w[i][j],mrf.w[i][j])
    return w_score


def compute_selfscore(mrf, edges_map, v_score_function, w_score_function):
    v_score = sum([v_score_function(vi,vi) for vi in mrf.v])
    w_score = compute_self_w_scores(mrf, edges_map, w_score_function)
    selfcomp = v_score+w_score
    return selfcomp
