import os
import csv
import numpy as np

from comutils import files_management as fm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_length_ins_file(input_file):
    """ returns MSA length (i.e. first sequence length without inserted lower case letters) """
    first_sequence = str(list(SeqIO.parse(str(input_file),'fasta'))[0].seq)
    return sum([1 for letter in first_sequence if not letter.islower()])


def count_insertions(input_file):
    """ returns delta_ins[nb sequences, msa length] where delta_ins[n,k] = length of gap at pos k for sequence n 
    @pc_insertions_tau is the pseudo-count coefficient (the higher tau, the more gap lengths are close to background probabilities for gap lengths)"""
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




def maximize_likelihood(delta_ins, L, maxit_infer_insertions=1e8, tol_infer_insertions=1e-6, learning_coeff_insertions=1e-3, freq_insert_min=1e-3, pc_insertions_tau=0, delta_n_max=100, **kwargs):
    """ maximizes log-likelihood to compute optimal gap open and gap extend penalties for each position
    returns a dictionary {"open":[list of gap open penalties], "extend":[list of gap extend penallties]}"""
    
    nseq = delta_ins.shape[0]
    expected_Nt = nseq*sum([delta_n*get_background_gap_probability(delta_n) for delta_n in range(delta_n_max)])
    expected_No = nseq*sum([get_background_gap_probability(delta_n) for delta_n in range(delta_n_max)])

    insertion_penalties = {"open":[0]*(L+1), "extend":[0]*(L+1)}


    for pos in range(L+1):

        eps=1e8
        lo=1.0
        le=1.0
        it=1

        No = sum(delta_ins[:,pos]>=1) # nb gap opens at pos (TODO vectorize) # >=1 ?? TODO réfléchir
        Nt = sum(delta_ins[:,pos]) # sum of lengths of gaps at pos (TODO vectorize)

        if (No==0):
            No=freq_insert_min*nseq
            Nt=freq_insert_min*nseq

        No = (1-pc_insertions_tau)*No+pc_insertions_tau*expected_No
        Nt = (1-pc_insertions_tau)*Nt+pc_insertions_tau*expected_Nt

    
        while ( (eps > tol_infer_insertions) and (it < maxit_infer_insertions) ):
            dLdlo = -No+nseq*(np.exp(-lo)/(1-np.exp(-le)+np.exp(-lo)))-2*lo
            lo = lo + learning_coeff_insertions * dLdlo
            dLdle = No-Nt+nseq*(np.exp(-lo-le)/((1-np.exp(-le))*(1-np.exp(-le)+np.exp(-lo))))-2*le
            le = le + learning_coeff_insertions * dLdle
            eps = max(abs(dLdle), abs(dLdlo))
            it+=1
       
        insertion_penalties["open"][pos] = lo
        insertion_penalties["extend"][pos] = le

    return insertion_penalties



def get_background_gap_probability(delta_n):
    # from Qian & Goldstein, 2001
    return 1.027e-2*np.exp(-delta_n/0.96)+3.031e-3*np.exp(-delta_n/3.13)+6.141e-4*np.exp(-delta_n/14.3)+2.090e-5*np.exp(-delta_n/81.7)



def infer_insertion_penalties(input_file, free_end_gaps=True, pc_insertions_tau=0, **kwargs):
    delta_ins = count_insertions(input_file)
    L = get_length_ins_file(input_file)
    insertion_penalties = maximize_likelihood(delta_ins, L, pc_insertions_tau=pc_insertions_tau, **kwargs)
    if free_end_gaps:
        for insertion_type in ['open', 'extend']:
            insertion_penalties[insertion_type][0] = 0
            insertion_penalties[insertion_type][-1] = 0
    return insertion_penalties



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



def infer_insertion_penalties_in_file(seed_a3m_file, output_file, free_end_gaps=True):
    """ calls DCAbuild function to infer insertion penalties from MSA @seed_a3m_file with insertions indicated as lower case letters wrt MSA of length @seed_length, output .tsv file in DCAbuild format @output_file"""
    fm.check_if_file_ok(seed_a3m_file)
    insertion_penalties = infer_insertion_penalties(seed_a3m_file, free_end_gaps=free_end_gaps)
    write_insertion_penalties_in_file(insertion_penalties, output_file)



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
