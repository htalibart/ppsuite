import os
from urllib.request import urlopen
import re
import json
import pandas as pd
from Bio import SeqIO, AlignIO, pairwise2
import Bio.PDB
from Bio.PDB.Polypeptide import PPBuilder
import ctypes
import pathlib
import shutil
import csv


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


def get_nb_sequences_in_fasta_file(fasta_file):
    records = list(SeqIO.parse(str(fasta_file), "fasta"))
    return len(records)


def create_fasta_file_with_less_sequences(aln_file, aln_1000, nb_sequences=1000, fileformat="fasta"):
    records = list(SeqIO.parse(str(aln_file),fileformat))
    with open(str(aln_1000), 'w') as f:
        SeqIO.write(records[:nb_sequences], f, fileformat)


def split_fasta(fasta_file, seq_folder):
    """ splits @fasta_file into multiple fasta files with 1 single sequence in folder @øeq_folder """
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
        shortest = pathlib.Path(files[0])
        if len(files)>1:
            print("more than 1 file ending with "+extension+", using the one with the shortest name : "+files[0])
        return shortest
    else:
        return None

def get_potts_model_file_from_folder(folder, mrf_type=None):
    if mrf_type is not None:
        p = get_file_from_folder_ending_with_extension(folder, "_"+mrf_type+".mrf")
        if p is not None:
            return p
    else:
        return get_file_from_folder_ending_with_extension(folder, ".mrf")

def get_sequence_file_from_folder(folder):
   return get_file_from_folder_ending_with_extension(folder, ".fasta")

def get_a3m_file_from_folder(folder):
   return get_file_from_folder_ending_with_extension(folder, ".a3m")

def get_pdb_file_from_folder(folder):
    pdb_file = get_file_from_folder_ending_with_extension(folder, ".pdb")
    if pdb_file is None:
        pdb_file = get_file_from_folder_ending_with_extension(folder, ".cif")
    return pdb_file


def write_readme(folder, **kwargs):
    p = folder/'README.txt'
    with p.open(mode='w') as f:
        json.dump(kwargs, f, default=str)

def copy(old_location, new_location):
    if not new_location.is_file():
        shutil.copy(old_location, new_location)

def get_format(seq_file):
    extension = seq_file.suffix[1:]
    if extension=="fa":
        return "fasta"
    else:
        return extension


def get_list_from_csv(csv_file):
    with open(csv_file, 'r') as f:
        csvreader = csv.reader(f)
        row = next(csvreader)
        l = [int(s) for s in row]
    return l

def write_list_to_csv(l, csv_file):
    with open(csv_file, 'w') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(l)

def fetch_pdb_file(pdb_id, outputfname):
    try:
        url = "https://files.rcsb.org/download/"+pdb_id+".pdb"
        pdbfile = urlopen(url)
        with open(outputfname+".pdb",'wb') as output:
            output.write(pdbfile.read())
        return outputfname+".pdb"
    except Exception as e:
        url = "https://files.rcsb.org/download/"+pdb_id+".cif"
        ciffile = urlopen(url)
        with open(outputfname+".cif", 'wb') as output:
            output.write(ciffile.read())
        return outputfname+".cif"

def get_pdb_chain(pdbid, pdb_file, chain_id='A'):
    pdbfile = str(pdb_file)
    if pdbfile.endswith(".pdb"):
        structure = Bio.PDB.PDBParser().get_structure(pdbid, pdbfile)
    elif pdbfile.endswith(".cif"):
        structure = Bio.PDB.MMCIFParser().get_structure(pdbid, pdbfile)
    else:
        raise Exception("Unknown PDB file format")
    model = structure[0]
    chain = model[chain_id]
    return chain

def get_sequence_from_pdb_chain(pdb_chain):
    ppb = PPBuilder()
    pdb_sequence = ppb.build_peptides(pdb_chain)[0].get_sequence()
    return pdb_sequence

def check_if_file_ok(f):
    if f is None:
        raise Exception(str(f)+" is not defined")
    elif not os.path.exists(str(f)):
        raise Exception("File not found :"+str(f))

def write_positions_to_csv(positions_dict, output_file):
    with open(str(output_file), 'w') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerow(list(positions_dict.keys()))
        for ind in range(len(positions_dict[list(positions_dict.keys())[0]])):
            row = [positions_dict[key][ind] for key in positions_dict]
            csvwriter.writerow(row)
