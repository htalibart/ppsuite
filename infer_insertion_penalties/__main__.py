import numpy as np
import ctypes
import os
import csv
import math
import pathlib
import tempfile

from comutils import files_management as fm
from comutils.tool_wrapper import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pkg_resources
INFER_INSERTIONS_CPP_LIBRARY_PATH = pkg_resources.resource_filename('infer_insertion_penalties', 'inference_lib.so')
INFERENCE_LIB = ctypes.CDLL(INFER_INSERTIONS_CPP_LIBRARY_PATH)



def get_length_ins_file(input_file):
    """ returns MSA length (i.e. first sequence length without inserted lower case letters) """
    first_sequence = str(list(SeqIO.parse(str(input_file),'fasta'))[0].seq)
    return sum([1 for letter in first_sequence if not letter.islower()])


#def count_insertions(input_file): # from DCAbuild
#    with open(input_file, 'r') as f:
#        records = list(SeqIO.parse(str(input_file), 'fasta'))
#        L = get_length_ins_file(input_file)
#        nseq = len(records)
#
#        ins = np.zeros((nseq,L+1))
#        
#        pointseq = 0
#
#        for rec in records:
#                s = str(rec.seq)
#                pointseq += 1
#                pos = L
#                for k in range(0,len(s)-1):
#                        if (s[len(s)-k-1].isupper()) or (s[len(s)-k-1] == '-'):
#                                go = True
#                                j = 1
#                                delta = 0
#                                while go:
#                                        if (len(s)-k-j <= 0):
#                                                go = False
#                                        elif  (s[len(s)-k-j-1].isupper()) or (s[len(s)-k-j-1] == '-'):
#                                                go = False
#                                        elif s[len(s)-k-j-1] == '.':
#                                                j += 1
#                                        elif s[len(s)-k-j-1].islower():
#                                                delta += 1
#                                                j += 1
#                                ins[pointseq-1,pos-1] = delta
#                                pos -= 1
#
#        return ins
#



def count_insertions(input_file):
    """ returns delta_ins[nb sequences, msa length] where delta_ins[n,k] = length of gap at pos k for sequence n""" 
    with open(input_file, 'r') as f:
        records = list(SeqIO.parse(str(input_file), 'fasta'))
        L = get_length_ins_file(input_file)
        nseq = len(records)

        delta_ins = np.zeros((nseq,L+1))

        for rec_index, rec in enumerate(records):

            sequence = str(rec.seq)
            in_gap=False

            index_in_msa=-1

            for letter in sequence:

                if letter.islower():
                    if not in_gap:
                        in_gap=True
                    delta_ins[rec_index,index_in_msa+1]+=1
                    
                else:
                    in_gap=False
                    index_in_msa+=1
            
        return delta_ins




def maximize_likelihood(delta_ins, length, output_insertions_file, maxit_infer_insertions=100000000, tol_infer_insertions=1e-6, learning_coeff_insertions=1e-3, freq_insert_min=1e-3, pc_insertions_tau=0, delta_n_max=100, **kwargs):

    c_double_p = ctypes.POINTER(ctypes.c_double) # pointer to double
    delta_ins_flat = np.ascontiguousarray(delta_ins.flatten()) # flatten array for ctypes
    c_delta_ins = delta_ins_flat.astype(np.float64).ctypes.data_as(c_double_p)


    nseq = delta_ins.shape[0]

    expected_No = nseq*sum([get_background_gap_probability(delta_n) for delta_n in range(delta_n_max)])
    expected_Nt = nseq*sum([delta_n*get_background_gap_probability(delta_n) for delta_n in range(delta_n_max)])


    #freq_insert_min=min(freq_insert_min,1/nseq)

    INFERENCE_LIB.call_infer_insertion_penalties_from_python.argtypes=[
            c_double_p, # delta_ins
            ctypes.c_int, # length
            ctypes.c_int, # nseq
            ctypes.c_int, # maxit
            ctypes.c_double, # tol
            ctypes.c_double, # learning_coeff
            ctypes.c_double, # freq_insert_min
            ctypes.c_double, # pc_insertions_tau
            ctypes.c_int, # delta_n_max
            ctypes.c_double, # expected No
            ctypes.c_double, # expected Nt
            ctypes.c_char_p # output_insertions_filename
            ]

    INFERENCE_LIB.call_infer_insertion_penalties_from_python(
            c_delta_ins,
            ctypes.c_int(length),
            ctypes.c_int(nseq),
            ctypes.c_int(maxit_infer_insertions),
            ctypes.c_double(tol_infer_insertions),
            ctypes.c_double(learning_coeff_insertions),
            ctypes.c_double(freq_insert_min),
            ctypes.c_double(pc_insertions_tau),
            ctypes.c_int(delta_n_max),
            ctypes.c_double(expected_No),
            ctypes.c_double(expected_Nt),
            ctypes.c_char_p(str(output_insertions_file).encode('utf-8'))
            )


def get_background_gap_probability(delta_n):
    # from Qian & Goldstein, 2001
    return 1.027e-2*math.exp(-delta_n/0.96)+3.031e-3*math.exp(-delta_n/3.13)+6.141e-4*math.exp(-delta_n/14.3)+2.090e-5*math.exp(-delta_n/81.7)



def infer_insertion_penalties_in_file(input_file, output_file, pc_insertions_tau=0, learning_coeff_insertions=1e-3, **kwargs):
    fm.check_if_file_ok(input_file)
    delta_ins = count_insertions(input_file)
    L = get_length_ins_file(input_file)
    maximize_likelihood(delta_ins, L, output_file, learning_coeff_insertions=learning_coeff_insertions, pc_insertions_tau=pc_insertions_tau, **kwargs)



def get_insertion_penalties_from_file(insertion_penalties_file):
    """ reads insertion penalties files in .tsv format into a dictionary {"open":[list of gap open penalties], "extend":[list of gap extend penalties]}"""
    fm.check_if_file_ok(insertion_penalties_file)
    insertion_penalties = {'open':[], 'extend':[]}

    # read penalties from file
    with open(str(insertion_penalties_file), 'r') as tsv_file:
        csv_reader = csv.reader(tsv_file, delimiter='\t')
        for row in csv_reader:
            insertion_penalties['open'].append(float(row[0]))
            insertion_penalties['extend'].append(float(row[1]))

    return insertion_penalties


def write_insertion_penalties_in_file(insertion_penalties, output_file):
    """ writes dictionary of insertion penalties {"open":[list of gap open penalties], "extend":[list of gap extend penalties]} in .tsv file in DCAbuild format """
    nb_inser = len(insertion_penalties['open'])
    with open(str(output_file), 'w') as tsv_file:
        csv_writer = csv.writer(tsv_file, delimiter='\t')
        for i in range(nb_inser):
            csv_writer.writerow([insertion_penalties['open'][i], insertion_penalties['extend'][i]])



def infer_insertion_penalties_in_dict(input_file, pc_insertions_tau=0, learning_coeff_insertions=1e-3, **kwargs):
    output_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names())) 
    infer_insertion_penalties_in_file(input_file, output_file, pc_insertions_tau=pc_insertions_tau, learning_coeff_insertions=learning_coeff_insertions, **kwargs)
    insertion_penalties = get_insertion_penalties_from_file(output_file)
    output_file.unlink()
    return insertion_penalties



def lower_case_trimmed_columns(aln_with_insertions, output_file, list_of_columns_not_trimmed, fileformat='fasta'):
    """ lowers letters of columns in @aln_with_insertions that are not in @list_of_columns_not_trimmed, writes result in @output_file """
    records = list(SeqIO.parse(str(aln_with_insertions),fileformat))
    output_records = []
    for record_index, record in enumerate(records):
        upper_index=-1
        new_record = record
        sequence_str = str(record.seq)
        new_sequence_str=''
        for letter in sequence_str:
            if letter.isupper() or letter=='-':
                upper_index+=1
                if upper_index not in list_of_columns_not_trimmed:
                    if letter!='-':
                        letter = letter.lower()
                    else:
                        letter=''
            new_sequence_str+=letter
        new_record.seq = Seq(new_sequence_str)
        output_records.append(new_record)
    SeqIO.write(output_records, str(output_file), 'fasta')


def ins_file_to_trimmed_ins_file(ins_file, output_file, trimal_gt, trimal_cons=0):
    reformat_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names())) 
    call_reformat(ins_file, reformat_file)
    trimmed_fasta = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names())) 
    colnum_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names())) 
    call_trimal(reformat_file, trimmed_fasta, trimal_gt, trimal_cons, colnumbering_file=colnum_file)
    list_of_columns_not_trimmed = fm.get_trimal_ncol(colnum_file)
    lower_case_trimmed_columns(ins_file, output_file, list_of_columns_not_trimmed, fileformat='fasta')
    reformat_file.unlink()
    trimmed_fasta.unlink()
    colnum_file.unlink()


