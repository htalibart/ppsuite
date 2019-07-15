import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def get_real_aligned_positions(aligned_positions, compotts_objects): # TODO test
    real_aligned_positions = {}
    c_names = ['pos_ref', 'pos_2']
    for k in range(2):
        real_aligned_positions[c_names[k]] = compotts_objects[k].get_real_positions(aligned_positions[c_names[k]])
    return real_aligned_positions


def get_seqs_aligned(aligned_positions, compotts_objects): #TODO test
    dict_pos = get_real_aligned_positions(aligned_positions, compotts_objects)
    c_names = ['pos_ref', 'pos_2']
    seqs_aligned = {ck:"" for ck in c_names}

    for k, ck in zip(range(2), c_names):
        seqs_aligned[ck]+='-'*(max(dict_pos['pos_ref'][0], dict_pos['pos_2'][0]))+compotts_objects[k].real_seq[dict_pos[ck][0]]

    old = [dict_pos['pos_ref'][0], dict_pos['pos_2'][0]]
    for pos in range(1,len(dict_pos['pos_ref'])):
        gap_length = max([dict_pos[ck][pos]-old[k] for k in range(2)])-1
        for k, ck in zip(range(2), c_names):
            seqs_aligned[ck]+='-'*gap_length+compotts_objects[k].real_seq[dict_pos[ck][pos]]
        old = [dict_pos['pos_ref'][pos],dict_pos['pos_2'][pos]]
    gap_length = max([len(compotts_objects[k].real_seq)-old[k] for k in range(2)])
    for ck in c_names:
        seqs_aligned[ck]+='-'*gap_length
    return [seqs_aligned['pos_ref'], seqs_aligned['pos_2']]



# TODO AlignIO ?
def get_seqs_aligned_in_fasta_file(aligned_positions, compotts_objects, output_file):
    seqs_aligned = get_seqs_aligned(aligned_positions, compotts_objects)
    seq_records = [SeqRecord(Seq(s, IUPAC.protein), id=o.name, description='') for s,o in zip(seqs_aligned, compotts_objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))



def get_seq_trimmed_for_ref(aligned_positions, ref_object, query_object): #TODO test
    aligned_positions_in_ref = aligned_positions['pos_ref']
    trimmed_seq=""
    for pos in range(ref_object.mrf.ncol):
        if pos in aligned_positions_in_ref:
            ind = aligned_positions_in_ref.index(pos)
            pos_2 = aligned_positions['pos_2'][ind]
            trimmed_seq+=query_object.get_real_letter_at_trimmed_pos(pos_2)
        else:
            trimmed_seq+='-'
    return trimmed_seq


def get_pos_aligned_at_pos(aligned_positions_dict, pos): #TODO test
    d = aligned_positions_dict
    if pos in d['pos_ref']:
        return d['pos_2'][d['pos_ref'].index(pos)]
    else:
        print(str(pos)+" is not aligned")
        return None



def get_aligned_v_scores(aligned_positions_dict, v_scores): #TODO test
    aligned_v_scores = np.zeros(len(aligned_positions_dict['pos_ref']))
    pos=0
    for i,k in zip(aligned_positions['pos_ref'], aligned_positions['pos_2']):
        aligned_v_scores[pos] = v_scores[i][k]
        pos+=1
    return aligned_v_scores
