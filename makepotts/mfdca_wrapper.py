from pydca.meanfield_dca import meanfield_dca
import numpy as np
from comutils.global_variables import ALPHABET
from makepotts import potts_model

INFINITY=100000000

RES_TO_INT_MFDCA = {
        'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
        'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10,
        'M': 11, 'N': 12, 'P': 13, 'Q': 14, 'R': 15,
        'S': 16, 'T': 17, 'V': 18, 'W':19, 'Y':20,
        '-':21
    }

INT_TO_RES_MFDCA = { RES_TO_INT_MFDCA[res]-1 : res for res in RES_TO_INT_MFDCA}


def get_ccmpredpy_index(mfdca_index):
    return ALPHABET.find(INT_TO_RES_MFDCA[mfdca_index])

def fields_list_to_v(fields_list):
    L = len(fields_list)
    v = np.zeros((L,21))
    for c in fields_list:
        for a in range(20):
            v[c[0],get_ccmpredpy_index(a)] = c[1][a]
    return v

def couplings_list_to_w(couplings_list, L):
    w = np.zeros((L,L,21,21))
    for c in couplings_list:
        for a in range(20):
            for b in range(20):
                w[c[0][0],c[0][1],get_ccmpredpy_index(a),get_ccmpredpy_index(b)] = c[1][a*20+b]
                w[c[0][1],c[0][0],get_ccmpredpy_index(b),get_ccmpredpy_index(a)] = w[c[0][0],c[0][1],get_ccmpredpy_index(a),get_ccmpredpy_index(b)]
    return w


def apply_zero_sum_gauge(v, w):
    zv = np.zeros_like(v)
    zw = np.zeros_like(w)
    L = v.shape[0]
    q = v.shape[1]
    average_v = np.mean(v, axis=1)
    average_w = np.mean(w, axis=(2,3))
    average_w_a = np.mean(w, axis=2)
    average_w_b = np.mean(w, axis=3)
    for i in range(L):
        for a in range(q):
            zv[i,a] = v[i,a]-average_v[i]+np.sum([average_w_b[i,j,a]-average_w[i,j] for j in range(L) if j!=i])
    for i in range(L-1):
        for j in range(i+1,L):
            for a in range(q):
                for b in range(q):
                    zw[i,j,a,b] = w[i,j,a,b]-average_w_b[i,j,a]-average_w_a[i,j,b]+average_w[i,j]
    return zv, zw



def infer_parameters(aln_file, mfdca_pseudocount=0.5, wt_cutoff=0.8, zero_sum_gauge=True, **kwargs):
    mfdca_inst = meanfield_dca.MeanFieldDCA(
        str(aln_file),
        'protein',
        pseudocount=mfdca_pseudocount,
        seqid = wt_cutoff,
    )
    fields_mf, couplings_mf = mfdca_inst.compute_params(num_site_pairs=INFINITY, linear_dist=0)
    L = len(fields_mf)
    v = fields_list_to_v(fields_mf)
    w = couplings_list_to_w(couplings_mf, L)
    if zero_sum_gauge:
        v, w = apply_zero_sum_gauge(v, w)
    return v, w

