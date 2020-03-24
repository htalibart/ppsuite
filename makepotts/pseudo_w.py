import numpy as np
from comutils.util import *
from makepotts.potts_model import *
from makepotts.rescaling import *
from makepotts.pair_substitution_matrices import P2P_PROBA, P2P_PROBA_CONTACT, P2P_PROBA_CONTACT_INTER, P2P_PROBA_CONTACT_INTER_NOT_WEIGHTED


def get_cond_proba(mat):
    background_prob = np.sum(mat, axis=1)[np.newaxis,:]
    cond = np.zeros_like(mat)
    for ab in range(len(mat)):
        for cd in range(len(mat)):
            cond[ab][cd] = mat[ab][cd]/background_prob[0][cd]
    return cond

def f_p2p(f, pair_matrix=P2P_PROBA):
    pb = np.sum(pair_matrix, axis=0)
    cond_prob = get_cond_proba(pair_matrix)
    q = 20
    f_pc = np.zeros_like(f)
    for a in range(q):
        for b in range(q):
            ab = a*q+b
            for c in range(q):
                for d in range(q):
                    cd = c*q+d
                    f_pc[:,:,a,b]+=cond_prob[ab,cd]*f[:,:,c,d]
    return f_pc


def f_with_pc(f, submat_tau, pair_matrix=P2P_PROBA):
    return (1-submat_tau)*np.copy(f)+submat_tau*f_p2p(f,pair_matrix=pair_matrix)


def get_direct_p(potts_model):
    # pij(a,b) = (1/Zij)exp(vi(a)+vj(b)+wij(a,b))
    p = np.zeros_like(potts_model.w)
    L = potts_model.ncol
    q = 20
    for i in range(L):
        for j in range(L):
            for a in range(q):
                for b in range(q):
                    p[i,j,a,b] = np.exp(potts_model.v[i][a]+potts_model.v[j][b]+potts_model.w[i][j][a][b])
            p[i,j,:,:] = p[i,j,:,:]/np.sum(p[i,j])
    return p


def p_to_w(p, potts_model):
    l = np.vectorize(almost_log)(p)
    w = np.zeros_like(p)
    L = potts_model.ncol
    q = 20
    for i in range(L):
        for j in range(L):
            for a in range(q):
                for b in range(q):
                    w[i,j,a,b] = l[i,j,a,b]-potts_model.v[i,a]-potts_model.v[j,b]-(1/(q*q))*np.sum(l[i][j])
    return w


def f_rescale_wij(wij, alpha_rescaling=10):
    resc_wij = np.vectorize(original_rescaling)(wij, alpha_rescaling=alpha_rescaling)
    if euclidean_norm(resc_wij)==0:
        return resc_wij
    else:
        return euclidean_norm(wij)/euclidean_norm(resc_wij)*resc_wij


def get_potts_model_with_pseudo_w(potts_model, w_submat_tau, rescale_wij=False, pair_matrix=P2P_PROBA):
    direct_p = get_direct_p(potts_model)
    direct_p_pc = f_with_pc(direct_p, w_submat_tau, pair_matrix=pair_matrix)
    w = p_to_w(direct_p_pc, potts_model)
    if not rescale_wij:
        return Potts_Model.from_parameters(potts_model.v, w)
    else:
        resc_w = np.zeros_like(w)
        L = potts_model.ncol
        for i in range(L):
            for j in range(L):
                resc_w[i][j] = f_rescale_wij(w[i][j])
        return Potts_Model.from_parameters(potts_model.v, resc_w)


def softmax_w(w):
    p = np.zeros_like(w)
    p[:,:,:,:] = np.exp(w[:,:,:,:])
    for i in range(len(p)):
        for j in range(len(p)):
            p[i,j,:,:] = p[i,j,:,:]/np.sum([p[i,j,:,:]])
    return p


def unsoftmax(p):
    w = np.zeros_like(p)
    l = np.vectorize(almost_log)(p)
    q = len(w[0][0])
    for i in range(len(w)):
        for j in range(len(w)):
            w[i,j,:,:] = l[i,j,:,:]-(1/(q*q))*np.sum(l[i,j,:,:])
    return w


def add_p2p_through_softmax(w, w_submat_tau=0.2):
    reduced_w = w[:,:,:20,:20]
    p = softmax_w(reduced_w)
    p_p2p = f_with_pc(p, submat_tau=w_submat_tau)
    extended_w = np.zeros_like(w)
    extended_w[:,:,:20,:20] = unsoftmax(p_p2p)
    return extended_w

def redistribute_probas(cond_matrix, alpha):
    reweighted_cond_matrix = np.zeros_like(cond_matrix)
    q = 20
    for ab in range(q*q):
        for cd in range(q*q):
            if ab!=cd:
                reweighted_cond_matrix[ab][cd] = cond_matrix[ab][cd]+(alpha/(q*q-1))*cond_matrix[cd][cd]
            else:
                reweighted_cond_matrix[ab][cd] = (1-alpha)*cond_matrix[ab][cd]
    return reweighted_cond_matrix

def reweight_wij(wij, cond_matrix):
    rwij = np.zeros_like(wij)
    q = 20
    for a in range(q):
        for b in range(q):
            rwij[a][b] = np.sum(wij*cond_matrix[a,b])
    if euclidean_norm(rwij)==0:
        return rwij
    else:
        return euclidean_norm(wij)/euclidean_norm(rwij)*rwij


def get_pseudo_wij(wij, reweighted_cond_matrix_4d, rescale, alpha_rescaling):
    rwij = reweight_wij(wij, reweighted_cond_matrix_4d)
    if rescale:
        return f_rescale_wij(rwij, alpha_rescaling)
    else:
        return rwij

def get_pseudo_w(w, reweighted_cond_matrix_4d, rescale, alpha_rescaling):
    pseudo_w = np.zeros_like(w)
    for i in range(len(w)):
        for j in range(len(w)):
            pseudo_w[i][j] = get_pseudo_wij(w[i][j], reweighted_cond_matrix_4d, rescale, alpha_rescaling)
    return pseudo_w

def add_pseudo_w_to_mrf(mrf, cond_matrix_2d=P2P_PROBA_CONTACT_INTER_NOT_WEIGHTED, alpha_probas=0.9, rescale=False, alpha_rescaling=30):
    reweighted_cond_matrix = matrix_to_4d(redistribute_probas(cond_matrix_2d, alpha_probas))
    return Potts_Model.from_parameters(mrf.v, get_pseudo_w(mrf.w, reweighted_cond_matrix, rescale, alpha_rescaling))
