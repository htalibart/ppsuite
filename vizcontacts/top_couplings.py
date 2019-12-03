#!/usr/bin/env python

""" Adapted from CCMpred https://github.com/soedinglab/CCMpred/blob/master/scripts/top_couplings.py """

import numpy as np


def main(f, output_file, reverse_order=False):
    mat = np.loadtxt(f)
    top = get_top_pairs(mat, reverse_order)
    of = open(output_file, "w")
    for i, j, coupling in zip(top[0], top[1], mat[top[1],top[0]]):
        of.write("{0},{1},{2}\n".format(i, j, coupling))
    of.close()


def get_top_pairs(mat, reverse_order):
    """Get the top-scoring contacts"""
    mat_masked = np.copy(mat)
    top = mat_masked.argsort(axis=None)
    if (not reverse_order):
        top = top[::-1]    
    top = (top % mat.shape[0]).astype(np.uint16), np.floor(top / mat.shape[0]).astype(np.uint16)
    return top

