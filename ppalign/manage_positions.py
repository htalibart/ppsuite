import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import comutils.files_management as fm
from comutils.util import *


def aln_dict_to_tuples_list(aln_dict):
    tuples_list = []
    for pos_aln in range(len(aln_dict["pos_ref"])):
        tuples_list.append((aln_dict["pos_ref"][pos_aln], aln_dict["pos_2"][pos_aln]))
    return tuples_list


def tuples_list_to_aln_dict(tuples_list):
    return {"pos_ref":[pair[0] for pair in tuples_list], "pos_2":[pair[1] for pair in tuples_list]}


def get_alignment_with_gaps(aligned_positions):
    tuples_list = aln_dict_to_tuples_list(aligned_positions)
    tuples_list_with_gaps = []
    previous_pair = (-1,-1)
    for pair in tuples_list:
        for pos_in_ref in list(range(previous_pair[0]+1,pair[0])):
            tuples_list_with_gaps.append((pos_in_ref,'-'))
        for pos_in_2 in list(range(previous_pair[1]+1,pair[1])):
            tuples_list_with_gaps.append(('-',pos_in_2))
        tuples_list_with_gaps.append(pair)
        previous_pair=pair
    return tuples_list_to_aln_dict(tuples_list_with_gaps)



def get_seq_positions_from_aln_dict(aligned_positions, objects):
    """ input : dict of lists of positions aligned by solver, output : dict of lists of positions in the original sequence """
    real_seq_positions = {}
    c_names = ['pos_ref', 'pos_2']
    for k in range(2):
        real_seq_positions[c_names[k]] = objects[k].get_seq_positions(aligned_positions[c_names[k]])
    real_seq_positions_without_none = {'pos_ref':[], 'pos_2':[]}
    for ind in range(len(real_seq_positions["pos_ref"])):
        if (real_seq_positions["pos_ref"][ind] is not None) and (real_seq_positions["pos_2"][ind] is not None):
            for name in c_names:
                real_seq_positions_without_none[name].append(real_seq_positions[name][ind])
    return real_seq_positions_without_none



def aligned_positions_to_aligned_sequences(seq_positions, sequences):
    c_names = ["pos_ref", "pos_2"]
    seqs_aligned = ["",""]
    for k in range(2):
        ck = c_names[k]
        for pos in seq_positions[ck]:
            if (pos=='-'):
                car='-'
            else:
                car=sequences[k][pos]
            seqs_aligned[k]+=car

    if (len(sequences[0])-seq_positions['pos_ref'][-1]+1>0):
        for pos_in_ref in range(seq_positions['pos_ref'][-1]+1,len(sequences[0])):
            seqs_aligned[0]+=sequences[0][pos_in_ref]
            seqs_aligned[1]+='-'
    if (len(sequences[1])-seq_positions['pos_2'][-1]+1>0):
        for pos_in_2 in range(seq_positions['pos_2'][-1]+1, len(sequences[1])):
            seqs_aligned[1]+=sequences[1][pos_in_2]
            seqs_aligned[0]+='-'
    return seqs_aligned


def get_seq_positions(aligned_positions, objects):
    seq_positions = {}
    seq_positions['pos_ref'] = objects[0].get_seq_positions(aligned_positions['pos_ref'])
    seq_positions['pos_2'] = objects[1].get_seq_positions(aligned_positions['pos_2'])
    return seq_positions


def remove_None_positions(positions_dict):
    new_positions_dict = {'pos_ref':[], 'pos_2':[]}
    for pos_aln in range(len(positions_dict['pos_ref'])):
        if (positions_dict['pos_ref'][pos_aln] is not None) and (positions_dict['pos_2'][pos_aln]) is not None:
            new_positions_dict['pos_ref'].append(positions_dict['pos_ref'][pos_aln])
            new_positions_dict['pos_2'].append(positions_dict['pos_2'][pos_aln])
    return new_positions_dict


def get_seqs_aligned(aligned_positions, objects):
    seq_positions_from_objects = get_seq_positions(aligned_positions, objects)
    seq_positions_without_None = remove_None_positions(seq_positions_from_objects)
    seq_positions = get_alignment_with_gaps(seq_positions_without_None)
    return aligned_positions_to_aligned_sequences(seq_positions, [obj.sequence for obj in objects]) 

            
def get_seqs_aligned_in_fasta_file(aligned_positions, objects, output_file):
    """ (positions aligned by solver + objects) -> sequences aligned -> in output_file """
    seqs_aligned = get_seqs_aligned(aligned_positions, objects)
    seq_records = [SeqRecord(Seq(s, IUPAC.protein), id=o.get_name(), description='') for s,o in zip(seqs_aligned, objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))



def aln_csv_to_aln_fasta(aln_file, objects, output_file):
    get_seqs_aligned_in_fasta_file(fm.get_aligned_positions_dict_from_ppalign_output_file(aln_file), objects, output_file)



def aln_sequences_csv_to_aln_fasta(aln_sequences_file, objects, output_file):
    aligned_positions = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_sequences_file)
    alignment_with_gaps = get_alignment_with_gaps(aligned_positions)
    aligned_sequences = aligned_positions_to_aligned_sequences(alignment_with_gaps, [obj.sequence for obj in objects])
    seq_records = [SeqRecord(Seq(s, IUPAC.protein), id=o.get_name(), description='') for s,o in zip(aligned_sequences, objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))

def aln_sequences_csv_to_aln_fasta_using_sequences_only(aln_sequences_file, sequence_files, output_file):
    aligned_positions = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_sequences_file)
    aligned_positions_without_None = remove_None_positions(aligned_positions)
    alignment_with_gaps = get_alignment_with_gaps(aligned_positions_without_None)
    records = [list(SeqIO.parse(seq_file, 'fasta'))[0] for seq_file in sequence_files]
    aligned_sequences = aligned_positions_to_aligned_sequences(alignment_with_gaps, [str(rec.seq) for rec in records])
    seq_records = [SeqRecord(Seq(s, IUPAC.protein), id=rec.id, description='') for s,rec in zip(aligned_sequences, records)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))



def get_pos_aligned_at_pos(aligned_positions_dict, pos):
    """ returns the position in object 2 aligned at position @pos in object 1"""
    d = aligned_positions_dict
    if pos in d['pos_ref']:
        return d['pos_2'][d['pos_ref'].index(pos)]
    else:
        print(str(pos)+" is not aligned")
        return None



def get_initial_positions(aligned_positions, mrf_pos_to_initial_pos_dict):
    initial_positions = {}
    for key in aligned_positions:
        initial_positions[key] = [mrf_pos_to_initial_pos_dict[key][pos] for pos in aligned_positions[key]]
    return initial_positions


def get_aln_sequences_from_aln_file(aln_file, objects, output_file):
    aligned_positions = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_file)
    sequence_positions = get_initial_positions(aligned_positions, {"pos_ref":objects[0].mrf_pos_to_seq_pos, "pos_2":objects[1].mrf_pos_to_seq_pos})
    fm.write_positions_to_csv(sequence_positions, output_file)

def get_mrf_pos_to_seq_pos(original_first_seq, seq, mrf_pos_to_aln_pos):
    seq_aln_pos = get_pos_first_seq_to_second_seq(original_first_seq, seq)
    mrf_pos_to_seq_pos = [seq_aln_pos[pos] for pos in mrf_pos_to_aln_pos]
    return mrf_pos_to_seq_pos



