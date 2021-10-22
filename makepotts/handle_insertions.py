import os
import csv
import numpy as np
import math

from comutils import files_management as fm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def maximize_likelihood(delta_ins, L, maxit_infer_insertions=1e8, tol_infer_insertions=1e-6, learning_coeff_insertions=1e-3, freq_insert_min=1e-3, pc_insertions_tau=0, delta_n_max=100, **kwargs): # deprecated
    """ maximizes log-likelihood to compute optimal gap open and gap extend penalties for each position
    returns a dictionary {"open":[list of gap open penalties], "extend":[list of gap extend penallties]}"""
    
    nseq = delta_ins.shape[0]
    expected_Nt = nseq*sum([delta_n*get_background_gap_probability(delta_n) for delta_n in range(delta_n_max)])
    expected_No = nseq*sum([get_background_gap_probability(delta_n) for delta_n in range(delta_n_max)])

    insertion_penalties = {"open":[0]*(L+1), "extend":[0]*(L+1)}


    le0 = None
    lo0 = None


    for pos in range(L+1):

        eps=1e8
        lo=1.0
        le=1.0
        it=1

        #No = sum(delta_ins[:,pos]>=1) # nb gap opens at pos (TODO vectorize) # >=1 ?? TODO réfléchir
        No = sum(delta_ins[:,pos]>0) # nb gap opens at pos (TODO vectorize) # >=1 ?? TODO réfléchir
        Nt = sum(delta_ins[:,pos]) # sum of lengths of gaps at pos (TODO vectorize)


        if (No==0):
            No=freq_insert_min*nseq
            no_gap=True
            #Nt=freq_insert_min*nseq
        else:
            no_gap=False


        #No = (1-pc_insertions_tau)*No+pc_insertions_tau*expected_No
        #Nt = (1-pc_insertions_tau)*Nt+pc_insertions_tau*expected_Nt

        fo=No/nseq
        ft=Nt/nseq


        if (fo<=freq_insert_min):
        #if (fo==0):
            no_gap=True
            fo==freq_insert_min
        else:
            no_gap=False


        if (no_gap) and (lo0 is not None) and (le0 is not None):
            lo = lo0
            le = le0

        else:
            while ( (eps > tol_infer_insertions) and (it < maxit_infer_insertions) ):
                emlo = math.exp(-lo)
                emle = math.exp(-le)
                emlole = math.exp(-lo-le)
                #dLdlo = -No+nseq*(math.exp(-lo)/(1-math.exp(-le)+math.exp(-lo)))-2*lo
                #dLdlo = -fo+(math.exp(-lo)/(1-math.exp(-le)+math.exp(-lo)))-2*lo/nseq
                dLdlo = -fo+(emlo/(1-emle+emlo))-2*lo/nseq
                #dLdlo = ( math.exp(-lo) * ((1.0-math.exp(-le))**(-1)) ) / (1.0 + math.exp(-lo)* ( (1.0 -math.exp(-le))**(-1) )) - fo - (2.0/nseq) *lo
                lo = lo + learning_coeff_insertions * dLdlo
                #dLdle = No-Nt+nseq*(math.exp(-lo-le)/((1-math.exp(-le))*(1-math.exp(-le)+math.exp(-lo))))-2*le
                #dLdle = fo-ft+(math.exp(-lo-le)/((1-math.exp(-le))*(1-math.exp(-le)+math.exp(-lo))))-2*le/nseq
                dLdle = fo-ft+emlole/((1-emle)*(1-emle+emlo))-2*le/nseq
                #dLdle = ( math.exp(-lo - le) *((1.0 - math.exp(-le))**(-2)) )/(1.0 + math.exp(-lo) * ( (1.0-math.exp(-le))**(-1) )) - ft + fo - (2.0/nseq)*le
                le = le + learning_coeff_insertions * dLdle
                eps = max(abs(dLdle), abs(dLdlo))
                it+=1
            print("pos=", pos, "fo=", fo, "ft=", ft, "it=", "{:.5e}".format(it), "eps=", eps, "lo=", lo, "le=", le, "dLdlo=",dLdlo, "dLdle=",dLdle)

        if (no_gap) and (lo0 is None) and (le0 is None):
            lo0 = lo
            le0 = le
           
        insertion_penalties["open"][pos] = lo
        insertion_penalties["extend"][pos] = le

    return insertion_penalties

