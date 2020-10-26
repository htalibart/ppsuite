from comutils.util import *
from comutils.global_variables import *
from makepotts.rescaling import *
import numpy as np
import math
from numpy import linalg as LA


def get_background_v0(v_rescaling_function_name="identity", rescale_removed_v0=False, **kwargs):
    f0 = np.array(AA_BACKGROUND_FREQUENCIES)
    logf0 = np.log(f0)
    vaa = logf0-(1/len(f0))*np.sum(logf0)
    v0 = np.append(vaa, [0])
    tiled_v0 = np.tile(v0, (1,1))
    if rescale_removed_v0:
        v0 = get_rescaled_parameters(tiled_v0, v_rescaling_function_name, **kwargs)
    return v0


def count_edges(edges_map):
    return sum(edges_map.flatten())

def get_w_threshold(mrf, w_percent):
    p = 100-int(w_percent)
    return np.percentile(mrf.get_w_norms(), p)


def get_edges_map(mrf, w_percent):
    w_threshold = get_w_threshold(mrf, w_percent)
    return 1*(mrf.get_w_norms()>w_threshold)

def get_vi_vk_score(vi, vk, remove_v0=False, offset_v=0, v_score_function=scalar_product, rescale_removed_v0=False, **kwargs):
    if remove_v0:
        v_bg = get_background_v0(rescale_removed_v0=rescale_removed_v0, **kwargs)
    else:
        v_bg=0
    return v_score_function(vi-v_bg,vk-v_bg)-offset_v


def compute_v_scores(mrf1, mrf2, v_score_function, offset_v, remove_v0, rescale_removed_v0, **kwargs):
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = get_vi_vk_score(mrf1.v[i], mrf2.v[k], remove_v0, offset_v, v_score_function, rescale_removed_v0, **kwargs)
    return v_scores


def get_wij_wkl_score(wij, wkl, w_score_function=scalar_product, **kwargs):
    return w_score_function(wij, wkl)


#def compute_w_scores(mrf1, mrf2, edges_map1, edges_map2, w_score_function, w_coeff=1, **kwargs):
#    """ symmetric matrix : w[i][j]=v[j+i*(i+1)/2])"""
#    len1 = int(mrf1.ncol*(mrf1.ncol+1)/2)
#    len2 = int(mrf2.ncol*(mrf2.ncol+1)/2)
#    print("computing w scores (LA="+str(mrf1.ncol)+", LB="+str(mrf2.ncol)+" : "+str(len1)+"x"+str(len2)+")")
#    w_scores = np.zeros((len1,len2), np.dtype('f4'))
#    for i in range(mrf1.ncol):
#        for j in range(i+1):
#            if edges_map1[i][j]:
#                for k in range(mrf2.ncol):
#                    for l in range(k+1):
#                        if edges_map2[k][l]:
#                            w_scores[j+int(i*(i+1)/2)][l+int(k*(k+1)/2)] = get_wij_wkl_score(mrf1.w[i][j], mrf2.w[k][l])
#    return w_coeff*w_scores
#
#
#def get_w_score_at(w_scores,i,k,j,l):
#    if (i==j):
#        return 0
#    elif (k==l):
#        return 0
#    else:
#        if (l<k):
#            if (j<i):
#                return w_scores[j+int(i*(i+1)/2)][l+int(k*(k+1)/2)]
#            else:
#                return get_w_score_at(w_scores,j,k,i,l)
#        else:
#            return get_w_score_at(w_scores,i,l,j,k)
#


def get_epsilon(epsilon_sim, selfscores):
    return (epsilon_sim/2)*sum(selfscores)


def compute_self_w_scores(mrf, edges_map, w_score_function, **kwargs):
    w_score = 0
    for i in range(mrf.ncol-1):
        for j in range(i+1,mrf.ncol):
            if edges_map[i][j]:
                w_score+=get_wij_wkl_score(mrf.w[i][j], mrf.w[i][j], w_score_function)
    return w_score


def compute_selfscore(mrf, edges_map, alpha_w=1, remove_v0=False, offset_v=0, use_v=True, use_w=True, v_score_function=scalar_product, w_score_function=scalar_product, rescale_removed_v0=False,**kwargs):
    if use_v:
        v_score = sum([get_vi_vk_score(vi, vi, remove_v0, offset_v, v_score_function=v_score_function, rescale_removed_v0=rescale_removed_v0, **kwargs) for vi in mrf.v])
    else:
        v_score = 0
    if use_w:
        w_score = compute_self_w_scores(mrf, edges_map, w_score_function)
    else:
        w_score = 0
    selfcomp = v_score+alpha_w*w_score
    return selfcomp


def get_v_scores_for_alignment(aligned_potts_models, aligned_positions_dict, remove_v0=False, offset_v=0, v_score_function=scalar_product, rescale_removed_v0=False, **kwargs):
    aligned_pos = [aligned_positions_dict["pos_ref"],aligned_positions_dict["pos_2"]]
    v_scores = np.array([get_vi_vk_score(aligned_potts_models[0].v[i],aligned_potts_models[1].v[k], remove_v0, offset_v, v_score_function, rescale_removed_v0=rescale_removed_v0, **kwargs) for i,k in zip(aligned_pos[0], aligned_pos[1])])
    return v_scores



def get_v_score_for_alignment(aligned_potts_models, aligned_positions_dict, remove_v0, offset_v, v_score_function=scalar_product, rescale_removed_v0=False, **kwargs):
    v_scores = get_v_scores_for_alignment(aligned_potts_models, aligned_positions_dict, remove_v0, offset_v, v_score_function, rescale_removed_v0=rescale_removed_v0, **kwargs)
    return np.sum(v_scores)


def get_w_scores_for_alignment(aligned_potts_models, dict_aligned_pos, w_score_function=scalar_product, **kwargs):
    aligned_pos = [dict_aligned_pos["pos_ref"],dict_aligned_pos["pos_2"]]
    L = len(aligned_pos[0])
    w_scores = np.zeros((L,L))
    for ind_i in range(L-1):
        for ind_j in range(ind_i+1,L):
            w_scores[ind_i,ind_j] = get_wij_wkl_score(aligned_potts_models[0].w[aligned_pos[0][ind_i],aligned_pos[0][ind_j]], aligned_potts_models[1].w[aligned_pos[1][ind_i],aligned_pos[1][ind_j]], w_score_function, **kwargs)
            w_scores[ind_j,ind_i]=w_scores[ind_i,ind_j]
    return w_scores


def get_w_score_for_alignment(aligned_potts_models, dict_aligned_pos, w_score_function=scalar_product, **kwargs):
   return 0.5*np.sum(get_w_scores_for_alignment(aligned_potts_models, dict_aligned_pos, w_score_function=w_score_function, **kwargs))


def get_total_gap_cost(ad, gap_open, **kwargs):
    gap_cost=0
    in_gap=True
    for pos_in_aln in range(1,len(ad["pos_ref"])):
        if not in_gap:
            if ((ad["pos_ref"][pos_in_aln]-ad["pos_ref"][pos_in_aln-1]>1) or (ad["pos_2"][pos_in_aln]-ad["pos_2"][pos_in_aln-1]>1)):
                in_gap=True
                gap_cost+=gap_open
        else:
            if ((ad["pos_ref"][pos_in_aln]-ad["pos_ref"][pos_in_aln-1]==1) and (ad["pos_2"][pos_in_aln]-ad["pos_2"][pos_in_aln-1]==1)):
                in_gap=False
    return gap_cost


def get_score_for_alignment(aligned_potts_models, aligned_positions_dict, alpha_w, **kwargs):
    return get_v_score_for_alignment(aligned_potts_models, aligned_positions_dict, **kwargs) + alpha_w * get_w_score_for_alignment(aligned_potts_models, aligned_positions_dict, **kwargs) - get_total_gap_cost(aligned_positions_dict, **kwargs)
