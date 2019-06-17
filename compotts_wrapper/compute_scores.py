from util import *
import numpy as np

# TODO vectoriser ou edges_map
# TODO décider d'un seuil

def count_edges(edges_map):
    nb=0
    for i in range(len(edges_map)):
        for j in range(len(edges_map)):
            nb+=edges_map[i][j]
    return nb


def get_edges_map(mrf, threshold):
    return 1*(mrf.get_w_norms()>threshold)


def compute_v_scores(mrf1, mrf2, v_score_function, **kwargs):
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = v_score_function(mrf1.v[i], mrf2.v[k])
    return v_scores


def compute_w_scores(mrf1, mrf2, edges_map1, edges_map2, w_score_function, **kwargs):
    print("computing w scores")
    w_scores = np.zeros((mrf1.ncol, mrf1.ncol, mrf2.ncol, mrf2.ncol))
    for i in range(mrf1.ncol-1):
        for j in range(i+1, mrf1.ncol):
            if edges_map1[i][j]:
                print(i,j)
                for k in range(mrf2.ncol-1):
                    for l in range(k+1, mrf2.ncol):
                        if edges_map2[k][l]: # TODO edges_map
                            w_scores[i][j][k][l] = w_score_function(mrf1.w[i][j], mrf2.w[k][l])
                            w_scores[j][i][k][l] = w_scores[i][j][k][l]
                            w_scores[i][j][l][k] = w_scores[i][j][k][l]
                            w_scores[j][i][l][k] = w_scores[i][j][k][l]
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
