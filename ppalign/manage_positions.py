import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import comutils.files_management as fm
from comutils.util import *

def get_alignment_with_gaps(aligned_positions, X='X'):
    """ input : dict of lists of aligned positions, output : alignment with gaps and "unknown areas" (symbole @X) """
    c_names = ["pos_ref", "pos_2"]
    aligned_positions_with_gaps = {ck:[] for ck in c_names}
    # positions 0 and before
    if aligned_positions["pos_ref"][0]==0:
        aligned_positions_with_gaps["pos_ref"]+=['-']*aligned_positions["pos_2"][0]+[0]
        aligned_positions_with_gaps["pos_2"]+=list(range(aligned_positions["pos_2"][0]+1))
    elif aligned_positions["pos_2"][0]==0:
        aligned_positions_with_gaps["pos_2"]+=['-']*aligned_positions["pos_ref"][0]+[0]
        aligned_positions_with_gaps["pos_ref"]+=list(range(aligned_positions["pos_ref"][0]+1))
    else:
        for ck in c_names:
            aligned_positions_with_gaps[ck]+=[X,aligned_positions[ck][0]]
    # positions [1:]
    for pos_aln in range(1,len(aligned_positions["pos_ref"])):
        diffs={ck : aligned_positions[ck][pos_aln]-aligned_positions[ck][pos_aln-1] for ck in c_names}
        if diffs["pos_ref"]==diffs["pos_2"]:
            for ck in c_names:
                aligned_positions_with_gaps[ck]+=list(range(aligned_positions[ck][pos_aln-1]+1,aligned_positions[ck][pos_aln]+1))
        elif diffs["pos_ref"]==1:
            aligned_positions_with_gaps["pos_ref"]+=['-']*(aligned_positions["pos_2"][pos_aln]-aligned_positions["pos_2"][pos_aln-1]-1)+[aligned_positions["pos_ref"][pos_aln]]
            aligned_positions_with_gaps["pos_2"]+=list(range(aligned_positions["pos_2"][pos_aln-1]+1,aligned_positions["pos_2"][pos_aln]+1))
        elif diffs["pos_2"]==1:
            aligned_positions_with_gaps["pos_2"]+=['-']*(aligned_positions["pos_ref"][pos_aln]-aligned_positions["pos_ref"][pos_aln-1]-1)+[aligned_positions["pos_2"][pos_aln]]
            aligned_positions_with_gaps["pos_ref"]+=list(range(aligned_positions["pos_ref"][pos_aln-1]+1,aligned_positions["pos_ref"][pos_aln]+1))
        else:
            for ck in c_names:
                aligned_positions_with_gaps[ck]+=[X]+[aligned_positions[ck][pos_aln]]
    return aligned_positions_with_gaps


def get_seq_positions(aligned_positions, objects):
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


def get_seqs_aligned(aligned_positions, objects, X='X'):
    """ input : dict of lists of aligned positions + corresponding objects, output : sequences aligned """
    c_names = ["pos_ref", "pos_2"]
    seq_positions = get_alignment_with_gaps(get_seq_positions(aligned_positions, objects), X=X)
    seqs_aligned = ["",""]
    for k in range(2):
        ck = c_names[k]
        for pos in seq_positions[ck]:
            if (pos=='-') or (pos=='X'):
                car=pos
            else:
                car=objects[k].sequence[pos]
            seqs_aligned[k]+=car
    return seqs_aligned

            
def get_seqs_aligned_in_fasta_file(aligned_positions, objects, output_file, X='X'):
    """ (positions aligned by solver + objects) -> sequences aligned -> in output_file """
    seqs_aligned = get_seqs_aligned(aligned_positions, objects, X=X)
    seq_records = [SeqRecord(Seq(s, IUPAC.protein), id=o.get_name(), description='') for s,o in zip(seqs_aligned, objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))


def get_seqs_aligned_in_fasta_file_only_aligned_positions(aligned_positions, objects, output_file):
    """ (positions aligned by solver + objects) -> sequences aligned (only at positions in the train msas) -> in output_file """
    msas = [list(SeqIO.parse(str(o.aln_train), "fasta")) for o in objects]
    records = []
    for k, c_name in zip(range(2), ['pos_ref', 'pos_2']):
        msa = msas[k]
        record = msa[0]
        new_record = record
        new_seq=""
        for pos in aligned_positions[c_name]:
            new_seq+=record[pos]
        new_record.seq=Seq(new_seq)
        records.append(new_record)
    with open(str(output_file), 'w') as f:
        SeqIO.write(records, f, "fasta")
    print("output can be found at "+str(output_file))




def get_pos_aligned_at_pos(aligned_positions_dict, pos):
    """ returns the position in object 2 aligned at position @pos in object 1"""
    d = aligned_positions_dict
    if pos in d['pos_ref']:
        return d['pos_2'][d['pos_ref'].index(pos)]
    else:
        print(str(pos)+" is not aligned")
        return None


def get_aligned_v_scores(aligned_positions_dict, v_scores):
    aligned_v_scores = np.zeros(len(aligned_positions_dict['pos_ref']))
    pos=0
    for i,k in zip(aligned_positions['pos_ref'], aligned_positions['pos_2']):
        aligned_v_scores[pos] = v_scores[i][k]
        pos+=1
    return aligned_v_scores


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



