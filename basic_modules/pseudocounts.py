import numpy as np
import ccmpred.substitution_matrices
import ccmpred.pseudocounts

SUBSTITUTION_MATRIX = ccmpred.substitution_matrices.BLOSUM62

def get_probas():
    return [np.sum(SUBSTITUTION_MATRIX[a]) for a in range(20)]
