import os
import re
import json
import pandas as pd
from Bio import SeqIO, AlignIO
import ctypes
import pathlib

def create_folder(name):
    p = pathlib.Path(name) 
    if not p.is_dir():
        p.mkdir()
    return p

def get_info_res_file_name(output_folder):
    return output_folder/"info.csv"

def get_aln_res_file_name(output_folder):
    return output_folder/"aln.csv"

def get_aligned_positions_dict_from_compotts_output_file(aln_res_file):
    df = pd.read_csv(aln_res_file)
    return df.to_dict('list') 

def get_infos_solver_dict_from_compotts_output_file(infos_res_file):
    df = pd.read_csv(infos_res_file)
    return df.loc[0].to_dict()


def create_seq_fasta(seq, fastaseq_file, seq_name="Billy"):
    """ crée un fichier fasta de nom @fastaseq_file avec la séquence @seq dedans, et retourne le nom de la séquence """
    with open(fastaseq_file, 'w') as of:
        of.write(">"+seq_name+"\n")
        of.write(seq+"\n")
    return seq_name


def get_sequence_by_name_in_fasta_file(seqname, fasta_file):
    record_dict = SeqIO.to_dict(SeqIO.parse(str(fasta_file), "fasta"))
    return record_dict[seqname]


def get_first_sequence_in_fasta_file(seq_file):
    return str(list(SeqIO.parse(str(seq_file), "fasta"))[0].seq)

def get_first_sequence_name(seq_file):
    return str(list(SeqIO.parse(str(seq_file), "fasta"))[0].id)


def get_first_sequence_clean_name(seq_file):
    actual_name = get_first_sequence_name(seq_file)
    return ''.join(e for e in actual_name if e.isalnum())


def get_nb_columns_in_alignment(aln_file):
    return len(get_first_sequence_in_fasta_file(aln_file))


def create_fasta_file_with_less_sequences(aln_file, aln_1000, nb_sequences=1000):
    AlignIO.write(AlignIO.read(str(aln_file), "fasta")[:nb_sequences], open(aln_1000, 'w'), "fasta")


def get_seq_names_from_seq_folder(seq_folder):
    return [str(f)[str(len(seq_folder)):-len(".fasta")] for f in seq_folder.glob('*.fasta')]

def separate_fasta(fasta_file, seq_folder):
    records = SeqIO.parse(open(str(fasta_file)), "fasta")
    create_folder(seq_folder)
    for record in records:
        with open(str(seq_folder)+str(record.id)+".fasta",'w') as f:
            SeqIO.write(record, f, "fasta")

def get_trimal_ncol(colnumbering_file):
    with colnumbering_file.open() as f:
        col_list = re.sub('[^0-9,.]', '', f.read()).split(',')
    return [int(s) for s in col_list]


def get_file_from_folder_ending_with_extension(folder, extension):
    files = [str(f) for f in pathlib.Path(folder).glob('*'+extension)]
    if len(files)>0:
        files.sort(key = len)
        shortest = files[0]
        if len(files)>1:
            print("more than 1 file ending with "+extension+", using the one with the shortest name : "+shortest)
        return shortest
    else:
        return None

def get_potts_model_file_from_folder(folder):
   return get_file_from_folder_ending_with_extension(folder, ".mrf")

def get_sequence_file_from_folder(folder):
   return get_file_from_folder_ending_with_extension(folder, ".fasta")

def get_a3m_file_from_folder(folder):
   return get_file_from_folder_ending_with_extension(folder, ".a3m")


def get_existing_training_set(folder, trimal_gt):
    trainf = get_file_from_folder_ending_with_extension(folder, "_training_set.fasta")
    if trainf is None:
        return get_file_from_folder_ending_with_extension(folder, "trim_"+str(int(trimal_gt*100))+".fasta")
    else:
        return trainf

def write_readme(folder, **kwargs):
    p = folder/'README.txt'
    with p.open(mode='w') as f:
        json.dump(kwargs, f, default=str)


