import numpy as np
from comutils.util import *
from makepotts.potts_model import *
from makepotts.rescaling import *
from makepotts.pair_substitution_matrices import P2P_PROBA

def f_p2p(f, pair_matrix=P2P_PROBA):
    pb = np.sum(pair_matrix, axis=0)
    cond_prob = pair_matrix / pb[np.newaxis, :]
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


def f_with_pc(f, submat_tau):
    return (1-submat_tau)*np.copy(f)+submat_tau*f_p2p(f)


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


def f_rescale_wij(wij):
    resc_wij = np.vectorize(original_rescaling)(wij, alpha_rescaling=10)
    if euclidean_norm(resc_wij)==0:
        return resc_wij
    else:
        return euclidean_norm(wij)/euclidean_norm(resc_wij)*resc_wij


def get_potts_model_with_pseudo_w(potts_model, w_submat_tau, rescale_wij=False):
    direct_p = get_direct_p(potts_model)
    direct_p_pc = f_with_pc(direct_p, w_submat_tau)
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
