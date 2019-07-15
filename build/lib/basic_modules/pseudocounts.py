import numpy as np
import ccmpred.substitution_matrices
import ccmpred.pseudocounts

SUBSTITUTION_MATRIX = ccmpred.substitution_matrices.BLOSUM62
# SUBSTITUTION_MATRIX[a][b] = p(a,b)


def get_cond_proba(a, knowing):
    p_knowing = np.sum(SUBSTITUTION_MATRIX[knowing])
    p_cond = SUBSTITUTION_MATRIX[a][knowing]/p_knowing
    return p_cond
