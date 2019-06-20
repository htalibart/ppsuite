import os
import re
import pandas as pd
from Bio import SeqIO, AlignIO
import ctypes

def get_compots_solver():
    return ctypes.CDLL(get_script_path()+"compotts_solver.so")

def create_folder(name):
    if name[-1]!='/':
        name+='/'
    if not os.path.isdir(name):
        os.mkdir(name)
    return name


def get_info_res_file_name(output_folder):
    return output_folder+"info.csv"

def get_aln_res_file_name(output_folder):
    return output_folder+"aln.csv"

def get_aligned_positions_dict_from_compotts_output_file(aln_res_file):
    df = pd.read_csv(aln_res_file)
    return df.to_dict('list') 

def get_infos_solver_dict_from_compotts_output_file(infos_res_file):
    df = pd.read_csv(infos_res_file)
    return df.loc[0].to_dict()


def create_seq_fasta(seq, fastaseq_file, seq_name="Billy"):
    """ crée un fichier fasta de nom @fastaseq_file avec la séquence @seq dedans, et retourne son nom """
    with open(fastaseq_file, 'w') as of:
        of.write(">"+seq_name+"\n")
        of.write(seq+"\n")
    return seq_name


def get_sequence_by_name_in_fasta_file(seqname, fasta_file):
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return record_dict[seqname]


def get_first_sequence_in_fasta_file(seq_file):
    return str(list(SeqIO.parse(seq_file, "fasta"))[0].seq)

def get_first_sequence_name(seq_file):
    return str(list(SeqIO.parse(seq_file, "fasta"))[0].id)


def get_name_from_first_sequence_name(seq_file):
    actual_name = get_first_sequence_name(seq_file)
    return ''.join(e for e in actual_name if e.isalnum())

def create_fasta_file_with_less_sequences(aln_file, aln_1000, nb_sequences=1000):
    AlignIO.write(AlignIO.read(aln_file, "fasta")[:nb_sequences], open(aln_1000, 'w'), "fasta")


def get_seq_names_from_seq_folder(seq_folder):
    seq_names = []
    for f in os.listdir(seq_folder):
        seq_names.append(f[:-len(".fasta")])
    return seq_names


def separate_fasta(fasta_file, seq_folder):
    records = SeqIO.parse(open(fasta_file), "fasta")
    create_folder(seq_folder)
    for record in records:
        with open(seq_folder+str(record.id)+".fasta",'w') as f:
            SeqIO.write(record, f, "fasta")

def get_trimal_ncol(colnumbering_file):
    with open(colnumbering_file, 'r') as f:
        col_list = re.sub('[^0-9,.]', '', f.read()).split(',')
    return [int(s) for s in col_list]


def get_script_path():
    return os.path.abspath(os.path.dirname(__file__))+'/' 
