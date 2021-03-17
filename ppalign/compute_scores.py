from comutils.util import *
from comutils.global_variables import *
from makepotts.rescaling import *
import numpy as np
import math
from numpy import linalg as LA


def get_background_v0(v_rescaling_function, rescale_removed_v0=False, **kwargs):
    """ returns background v0 vector field, rescaled with @v_rescaling_function if @rescaled_removed_v0 """
    f0 = np.array(AA_BACKGROUND_FREQUENCIES)
    logf0 = np.log(f0)
    vaa = logf0-(1/len(f0))*np.sum(logf0)
    v0 = np.append(vaa, [0])
    tiled_v0 = np.tile(v0, (1,1))
    if rescale_removed_v0:
        v0 = get_rescaled_parameters(tiled_v0, v_rescaling_function, **kwargs)
    return v0

def get_w_threshold(mrf, w_percent):
    """ returns the ||wij|| min considered by PPalign when considering the @w_percent % strongest couplings"""
    p = 100-int(w_percent)
    return np.percentile(mrf.get_w_norms(), p)


def get_edges_map(mrf, w_percent):
    """ returns a 2D array where array[i][j] if edge (i,j) should be considered by PPalign and 0 otherwise when considering the first @w_percent % strongest couplings """
    w_threshold = get_w_threshold(mrf, w_percent)
    return 1*(mrf.get_w_norms()>w_threshold)


def get_vi_vk_score(vi, vk, remove_v0=False, offset_v=0, v_score_function=scalar_product, rescale_removed_v0=False, **kwargs):
    """ returns the similarity score of field vectors @vi and @vk """
    if remove_v0:
        v_bg = get_background_v0(rescale_removed_v0=rescale_removed_v0, **kwargs)
    else:
        v_bg=0
    return v_score_function(vi-v_bg,vk-v_bg)-offset_v


def compute_v_scores(mrf1, mrf2, v_score_function, offset_v, remove_v0, rescale_removed_v0, **kwargs):
    """ returns an np.array @v_scores where v_scores[i][k] is the similarity score between field vi in Potts model @mrf1 and field vk in Potts model @mrf2 """
    v_scores = np.zeros((mrf1.ncol, mrf2.ncol))
    for i in range(mrf1.ncol):
        for k in range(mrf2.ncol):
            v_scores[i][k] = get_vi_vk_score(mrf1.v[i], mrf2.v[k], remove_v0, offset_v, v_score_function, rescale_removed_v0, **kwargs)
    return v_scores


def get_wij_wkl_score(wij, wkl, w_score_function=scalar_product, **kwargs):
    """ returns the similarity score of couplint matrices @wij and @wkl """
    return w_score_function(wij, wkl)


def get_epsilon(epsilon_sim, selfscores):
    """ gives epsilon for s(A,B) so that 2*s(A,B)/(s(A,A)+s(B,B)) < @epsilon_sim (where s(A,A) and s(B,B) are @selfscores)"""
    return (epsilon_sim/2)*sum(selfscores)


def compute_self_w_scores(mrf, edges_map, w_score_function, **kwargs):
    """ @w_score[i,j] is score for aligning w[i,j] with w[i,j] in Potts model @mrf for each considered edge in @edges_map """
    w_score = 0
    for i in range(mrf.ncol-1):
        for j in range(i+1,mrf.ncol):
            if edges_map[i][j]:
                w_score+=get_wij_wkl_score(mrf.w[i][j], mrf.w[i][j], w_score_function)
    return w_score


def compute_selfscore(mrf, edges_map, alpha_w=1, remove_v0=False, offset_v=0, use_v=True, use_w=True, v_score_function=scalar_product, w_score_function=scalar_product, rescale_removed_v0=False,**kwargs):
    """ score for aligning @mrf with itself """
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
    """ @v_scores[i] is the score for the alignment of field vectors at positions in the two @aligned_potts_models that are aligned at position i in the PPalign alignment given by @aligned_positions_dict """
    aligned_pos = [aligned_positions_dict["pos_ref"],aligned_positions_dict["pos_2"]]
    v_scores = np.array([get_vi_vk_score(aligned_potts_models[0].v[i],aligned_potts_models[1].v[k], remove_v0, offset_v, v_score_function, rescale_removed_v0=rescale_removed_v0, **kwargs) for i,k in zip(aligned_pos[0], aligned_pos[1])])
    return v_scores


def get_v_score_for_alignment(aligned_potts_models, aligned_positions_dict, remove_v0, offset_v, v_score_function=scalar_product, rescale_removed_v0=False, **kwargs):
    """ total v score """
    v_scores = get_v_scores_for_alignment(aligned_potts_models, aligned_positions_dict, remove_v0, offset_v, v_score_function, rescale_removed_v0=rescale_removed_v0, **kwargs)
    return np.sum(v_scores)


def get_w_scores_for_alignment(aligned_potts_models, dict_aligned_pos, w_score_function=scalar_product, **kwargs):
    """ @w_scores[i,j] is the score for the alignment of couplings from the two @aligned_potts_models that are aligned at position i and j in the PPalign alignment given by @aligned_positions_dict """
    aligned_pos = [dict_aligned_pos["pos_ref"],dict_aligned_pos["pos_2"]]
    L = len(aligned_pos[0])
    w_scores = np.zeros((L,L))
    for ind_i in range(L-1):
        for ind_j in range(ind_i+1,L):
            w_scores[ind_i,ind_j] = get_wij_wkl_score(aligned_potts_models[0].w[aligned_pos[0][ind_i],aligned_pos[0][ind_j]], aligned_potts_models[1].w[aligned_pos[1][ind_i],aligned_pos[1][ind_j]], w_score_function, **kwargs)
            w_scores[ind_j,ind_i]=w_scores[ind_i,ind_j]
    return w_scores


def get_w_score_for_alignment(aligned_potts_models, dict_aligned_pos, w_score_function=scalar_product, **kwargs):
    """ total w score (alpha_w not taken into account here """
    return 0.5*np.sum(get_w_scores_for_alignment(aligned_potts_models, dict_aligned_pos, w_score_function=w_score_function, **kwargs))


def get_total_gap_cost(ad, gap_open, sequence_lengths, **kwargs):
    """ inputs:
        @ad: aligned positions dict outputted by PPalign
        @sequence_lengths: list of lengths for each sequence
        outputs total gap cost for the alignment """
    gap_cost=0
    for counter, pos_type in enumerate(["pos_ref", "pos_2"]):
        in_gap=False
        previous_aln_pos=-1
        for pos_in_aln in range(len(ad[pos_type])):
            if not in_gap:
                if ad[pos_type][pos_in_aln]-previous_aln_pos>1:
                    in_gap=True
                    gap_cost+=gap_open
            else:
                if ad[pos_type][pos_in_aln]-previous_aln_pos==1:
                    in_gap=False
            previous_aln_pos = ad[pos_type][pos_in_aln]
        if ad[pos_type][-1]<(sequence_lengths[counter]-1):
            gap_cost+=gap_open
    return gap_cost


def get_score_for_alignment(aligned_potts_models, aligned_positions_dict, alpha_w, **kwargs):
    """ total score for PPalign alignment """
    return get_v_score_for_alignment(aligned_potts_models, aligned_positions_dict, **kwargs) + alpha_w * get_w_score_for_alignment(aligned_potts_models, aligned_positions_dict, **kwargs) - get_total_gap_cost(aligned_positions_dict, sequence_lengths=[mrf.ncol for mrf in aligned_potts_models], **kwargs)
