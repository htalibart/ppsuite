import csv
import numpy as np

import pkg_resources

def get_proba_mat_from_csv(csv_file):
    mat = np.zeros((400,400))
    with open(csv_file) as f:
        csvreader = csv.reader(f)
        i = 0
        for row in csvreader:
            for j in range(len(row)):
                mat[i][j] = float(row[j])
            i+=1
    return mat


P2P_PROBA = get_proba_mat_from_csv(pkg_resources.resource_filename(__name__,'P2Pmat14Weighted_prob.csv'))
P2P_PROBA_CONTACT = get_proba_mat_from_csv(pkg_resources.resource_filename(__name__,'P2Pmat14WeightedInterIntra=_prob_cond.csv'))
P2P_PROBA_CONTACT_INTER = get_proba_mat_from_csv(pkg_resources.resource_filename(__name__,'P2Pmat14WeightedInter=_prob_cond.csv'))
P2P_PROBA_CONTACT_INTER_NOT_WEIGHTED = get_proba_mat_from_csv(pkg_resources.resource_filename(__name__,'P2Pmat14Inter=_prob_cond.csv'))
