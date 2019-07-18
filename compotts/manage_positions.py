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


def get_seqs_aligned(aligned_positions, compotts_objects):
    real_aligned_positions = get_real_aligned_positions(aligned_positions, compotts_objects)
    c_names = ["pos_ref", "pos_2"]
    prec_pos = {ck : 0 for ck in c_names}
    seqs_aligned = {ck : "" for ck in c_names}
    seqs = {"pos_ref":compotts_objects[0].real_seq, "pos_2":compotts_objects[1].real_seq}
    for pos_aln in range(len(aligned_positions["pos_ref"])):
        pos = {}
        for ck in c_names:
            pos[ck] = real_aligned_positions[ck][pos_aln]
        diffs = {ck:pos[ck]-prec_pos[ck] for ck in c_names}
        if (diffs["pos_ref"]==diffs["pos_2"]):
            for ck in c_names:
                seqs_aligned[ck]+=seqs[ck][prec_pos[ck]+1:pos[ck]+1]
        else:
            if (diffs["pos_ref"]==0):
                seqs_aligned["pos_ref"]+="-"*diffs["pos_2"]
                seqs_aligned["pos_2"]+=seqs["pos_2"][prec_pos["pos_2"]+1:pos["pos_2"]+1]
            elif (diffs["pos_2"]==0):
                seqs_aligned["pos_2"]+="-"*diffs["pos_ref"]
                seqs_aligned["pos_ref"]+=seqs["pos_ref"][prec_pos["pos_ref"]+1:pos["pos_ref"]+1]
            else:
                for ck in c_names:
                    seqs_aligned[ck]+='X'+seqs[ck][pos[ck]]
        prec_pos = {ck : pos[ck] for ck in c_names}
    pos = {ck : len(seqs[ck])-1 for ck in c_names}
    for ck in c_names:
            pos[ck] = real_aligned_positions[ck][pos_aln]
    if (pos["pos_ref"]-prec_pos["pos_ref"])==(pos["pos_2"]-prec_pos["pos_2"]):
        for ck in c_names:
            seqs_aligned[ck]+=seqs[ck][prec_pos[ck]+1:pos[ck]+1]
    else:
        for ck in c_names:
            seqs_aligned[ck]+='X'
    return [seqs_aligned[ck] for ck in c_names]

            

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
