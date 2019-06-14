from util import *
import numpy as np

# TODO vectoriser


def compute_v_scores(mrf1, mrf2, v_score_function, **kwargs):
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = v_score_function(mrf1.v[i], mrf2.v[k])
    return v_scores


def compute_w_scores(mrf1, mrf2, w_score_function, **kwargs):
    w_scores = np.zeros((mrf1.ncol, mrf1.ncol, mrf2.ncol, mrf2.ncol))
    for i in range(mrf1.ncol-1):
        for j in range(i+1, mrf1.ncol):
            for k in range(mrf2.ncol-1):
                for l in range(k+1, mrf2.ncol):
                    print(i, j, k, l)
                    w_scores[i][j][k][l] = w_score_function(mrf1.w[i][j], mrf2.w[k][l])
                    w_scores[j][i][k][l] = w_scores[i][j][k][l]
                    w_scores[i][j][l][k] = w_scores[i][j][k][l]
                    w_scores[j][i][l][k] = w_scores[i][j][k][l]
    return w_scores

