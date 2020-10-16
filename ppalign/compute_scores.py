from comutils.util import *
import numpy as np
import math
from numpy import linalg as LA

def count_edges(edges_map):
    return sum(edges_map.flatten())


def get_w_threshold(mrf, w_percent):
    p = 100-int(w_percent)
    return np.percentile(mrf.get_w_norms(), p)


def get_edges_map(mrf, w_percent):
    w_threshold = get_w_threshold(mrf, w_percent)
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



def get_epsilon(epsilon_sim, selfscores):
    return (epsilon_sim/2)*sum(selfscores)


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
