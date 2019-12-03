import tempfile
import pathlib
import numpy as np

from collections import OrderedDict

from comfeature.comfeature import *
from vizcontacts import top_couplings

def get_contact_scores_for_aln_train(comfeature):

    # Get contact scores in a csv file
    mat_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
    apc_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
    ccmpredpy_call = "ccmpred "+str(comfeature.aln_train)+" -i "+str(comfeature.potts_model_file)+" --do-not-optimize -m "+str(mat_file)+" --apc "+str(apc_file)
    subprocess.Popen(ccmpredpy_call, shell=True).wait()
    output_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
    top_couplings.main(apc_file, output_file)

    # ordered dictionary {frozenset(c1,c2) : score}
    contact_scores = OrderedDict()
    with open(output_file, 'r') as f:
        csvreader = csv.reader(f)
        for row in csvreader:
            contact_scores[frozenset((int(row[0]),int(row[1])))] = float(row[2])

    interesting_contact_scores = OrderedDict()
    for c in contact_scores:
        if len(c)==2:
            interesting_contact_scores[tuple(c)] = contact_scores[c]

    return interesting_contact_scores


def get_contact_scores_for_sequence(comfeature):
    aln_scores = get_contact_scores_for_aln_train(comfeature)
    seq_contact_scores = OrderedDict()
    for c in aln_scores:
        c_seq = (comfeature.mrf_pos_to_seq_pos[c[0]], comfeature.mrf_pos_to_seq_pos[c[1]])
        seq_contact_scores[c_seq] = aln_scores[c]
    return seq_contact_scores


def translate_dict_to_pdb_pos(couplings_dict, pdb_chain, real_sequence):
    pdb_sequence = fm.get_sequence_from_pdb_chain(pdb_chain) 
    d = fm.get_pos_dict_first_seq_to_second_seq(real_sequence, pdb_sequence)
    pdb_couplings_dict = {}
    for c in couplings_dict:
        new_c = (d[c[0]], d[c[1]])
        if not None in new_c:
            pdb_couplings_dict[new_c] = couplings_dict[c]
    return pdb_couplings_dict


def is_true_contact(pdb_sequence_coupling, pdb_chain, contact_distance=8):
    return aa_distance(pdb_sequence_coupling[0], pdb_sequence_coupling[1], pdb_chain) <= contact_distance


def aa_distance(pos1, pos2, pdb_chain):
    r1 = pdb_chain[pos1]
    r2 = pdb_chain[pos2]
    diff_vector = r1['CA'].coord - r2['CA'].coord
    return np.sqrt(np.sum(diff_vector*diff_vector))


def get_colored_true_false_dicts(couplings_dict, pdb_chain, real_sequence, colors={True:'blue', False:'red'}, contact_distance=8):
    pdb_d = translate_dict_to_pdb_pos(couplings_dict, pdb_chain, real_sequence)
    tf_d = {colors[c]:{} for c in colors}
    for c in pdb_d:
        tf_d[colors[is_true_contact(c, pdb_chain, contact_distance=contact_distance)]][c] = pdb_d[c]
    return tf_d

def remove_couplings_too_close(couplings_dict, coupling_sep_min):
    ok_dict = {}
    for c in couplings_dict:
        if None not in c:
            if abs(c[0]-c[1])>=coupling_sep_min:
                ok_dict[c] = couplings_dict[c]
    return ok_dict
