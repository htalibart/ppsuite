import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def get_alignment_with_gaps(aligned_positions):
    """ input : dict of lists of aligned positions, output : alignment with gaps and "unknown areas" """
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
            aligned_positions_with_gaps[ck]+=['X',aligned_positions[ck][0]]
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
                aligned_positions_with_gaps[ck]+=['X']+[aligned_positions[ck][pos_aln]]
    print(aligned_positions_with_gaps)
    return aligned_positions_with_gaps


def get_seq_positions(aligned_positions, compotts_objects):
    """ input : dict of lists of positions aligned by solver, output : dict of lists of positions in the original sequence """
    real_seq_positions = {}
    c_names = ['pos_ref', 'pos_2']
    for k in range(2):
        real_seq_positions[c_names[k]] = compotts_objects[k].get_seq_positions(aligned_positions[c_names[k]])
    return real_seq_positions


def get_seqs_aligned(aligned_positions, compotts_objects):
    """ input : dict of lists of aligned positions + corresponding compotts objects, output : sequences aligned """
    c_names = ["pos_ref", "pos_2"]
    seq_positions = get_alignment_with_gaps(get_seq_positions(aligned_positions, compotts_objects))
    seqs_aligned = ["",""]
    for k in range(2):
        ck = c_names[k]
        for pos in seq_positions[ck]:
            if (pos=='-') or (pos=='X'):
                car=pos
            else:
                car=compotts_objects[k].sequence[pos]
            seqs_aligned[k]+=car
    return seqs_aligned

            

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


def get_real_pos_list(real_seq, other_seq):
    """ retourne une liste où other_to_real_dict[k] est la position dans la vraie séquence @real_seq correspondant à la position k de @other_seq """
    alignments = pairwise2.align.globalxx(other_seq, real_seq)
    aln = alignments[0][0]
    print(aln)
    pos_in_other_seq = 0
    pos_in_real_seq = 0
    other_to_real_dict = []
    for i in range(len(aln)):
        if (aln[i]!='-'):
            other_to_real_dict.append(pos_in_real_seq)
            pos_in_other_seq+=1
        pos_in_real_seq+=1
    return other_to_real_dict



def get_aligned_v_scores(aligned_positions_dict, v_scores): #TODO test
    aligned_v_scores = np.zeros(len(aligned_positions_dict['pos_ref']))
    pos=0
    for i,k in zip(aligned_positions['pos_ref'], aligned_positions['pos_2']):
        aligned_v_scores[pos] = v_scores[i][k]
        pos+=1
    return aligned_v_scores
