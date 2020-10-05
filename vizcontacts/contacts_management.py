import tempfile
import pathlib
import numpy as np
from collections import OrderedDict
from itertools import islice
from kneebow.rotor import Rotor
import matplotlib.pyplot as plt
from makepotts.potts_object import *
from comutils.util import *
from vizcontacts import top_couplings

def get_contact_scores_for_aln_train(potts_object):

    # Get contact scores in a csv file
    mat_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
    apc_file = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
    ccmpredpy_call = "ccmpred "+str(potts_object.aln_train)+" -i "+str(potts_object.potts_model_file)+" --do-not-optimize -m "+str(mat_file)+" --apc "+str(apc_file)
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


def get_contact_scores_for_sequence(potts_object):
    aln_scores = get_contact_scores_for_aln_train(potts_object)
    seq_contact_scores = OrderedDict()
    for c in aln_scores:
        c_seq = (potts_object.mrf_pos_to_seq_pos[c[0]], potts_object.mrf_pos_to_seq_pos[c[1]])
        seq_contact_scores[c_seq] = aln_scores[c]
    return seq_contact_scores


def get_pdb_offset(pdb_file, chain_id, real_sequence):
    pdb_chain = fm.get_pdb_chain(pdb_file, chain_id)
    r = next(pdb_chain.get_residues())
    offset = r.get_full_id()[3][1]-1
    return offset


def get_real_pos_to_pdb_pos(pdb_file, chain_id, real_sequence):
    pdb_sequence = fm.get_sequence_from_pdb_file(pdb_file, chain_id) 
    d = get_pos_first_seq_to_second_seq(real_sequence, pdb_sequence) # d[pos_in_real_seq] = pos_in_pdb_seq
    pdb_offset = get_pdb_offset(pdb_file, chain_id, real_sequence)
    rtpdb = []
    for pos in range(len(real_sequence)):
        if d[pos] is None:
            rtpdb.append(None)
        else:
            rtpdb.append(d[pos]+1+pdb_offset)
    return rtpdb


def translate_dict_to_pdb_pos(couplings_dict, pdb_file, chain_id, real_sequence):
    d = get_real_pos_to_pdb_pos(pdb_file, chain_id, real_sequence)
    pdb_couplings_dict = OrderedDict()
    for c in couplings_dict:
        if not None in c:
            new_c = (d[c[0]], d[c[1]])
            if not None in new_c:
                pdb_couplings_dict[new_c] = couplings_dict[c]
    return pdb_couplings_dict


def is_true_contact(pdb_sequence_coupling, pdb_file, chain_id, contact_distance=8):
    pdb_chain = fm.get_pdb_chain(pdb_file, chain_id)
    return aa_distance(pdb_sequence_coupling[0], pdb_sequence_coupling[1], pdb_chain) <= contact_distance


def aa_distance(pos1, pos2, pdb_chain):
    r1 = pdb_chain[pos1]
    r2 = pdb_chain[pos2]
    diff_vector = r1['CA'].coord - r2['CA'].coord
    return np.sqrt(np.sum(diff_vector*diff_vector))


def get_colored_true_false_dicts(couplings_dict, pdb_file, chain_id, real_sequence, colors={True:'blue', False:'red'}, contact_distance=8):
    d = get_real_pos_to_pdb_pos(pdb_file, chain_id, real_sequence)
    tf_d = {colors[val]:OrderedDict() for val in colors}
    for c in couplings_dict:
        pdb_c = (d[c[0]], d[c[1]])
        if not None in pdb_c:
            tf_d[colors[is_true_contact(pdb_c, pdb_file, chain_id, contact_distance=contact_distance)]][c] = couplings_dict[c]
    return tf_d

def remove_couplings_too_close(couplings_dict, coupling_sep_min):
    ok_dict = OrderedDict()
    for c in couplings_dict:
        if None not in c:
            if abs(c[0]-c[1])>coupling_sep_min:
                ok_dict[c] = couplings_dict[c]
    return ok_dict

def get_smaller_dict(couplings_dict, nb_couplings):
    new_dict = OrderedDict()
    for c in couplings_dict:
        if len(new_dict)<nb_couplings:
            new_dict[c] = couplings_dict[c]
    return new_dict

def get_elbow_index(couplings_dict, plot_elbow=False):
    y = list(couplings_dict.values())
    x = list(range(len(y)))

    if plot_elbow:
        plt.figure()
        plt.plot(x,y)
        plt.show()

    y.reverse()
    data = np.array([[xi,yi] for xi,yi in zip(x,y)])
    rotor = Rotor()
    rotor.fit_rotate(data)
    elbow_idx = rotor.get_elbow_index()

    return len(y)-elbow_idx

def get_cutoff_smaller_than(couplings_dict, score_cutoff):
    y = list(couplings_dict.values())
    if y[0]<score_cutoff:
        #return None
        return 0
    else:
        ind=0
        while (y[ind]>=score_cutoff):
            ind+=1
        return ind

def get_exclus_overlaps(couplings_dict_, tops):
    """ returns a list of 3 dicts : [couplings_dict[0]\couplings_dict[1], couplings_dict[0] inter couplings_dict[1], couplings_dict[1]\couplings_dict[0]. Overlap -> mean score """
    if len(couplings_dict_)==1:
        return [OrderedDict(islice(couplings_dict_[0].items(), 0, tops[0]))]
    else:
        couplings_dict = [OrderedDict(islice(couplings_dict_[k].items(), 0, int(tops[k]))) for k in range(2)]
        exclus = []
        overlaps = {}
        for k in range(2):
            exclus.append(OrderedDict())
            for c in couplings_dict[k]:
                if not ( (c in couplings_dict[(k+1)%2]) or ((c[1],c[0]) in couplings_dict[(k+1)%2]) ):
                    exclus[k][c] = couplings_dict[k][c]
                elif c in couplings_dict[(k+1)%2]:
                    overlaps[frozenset(c)] =  (couplings_dict[k][c]+couplings_dict[(k+1)%2][c])/2
                elif (c[1],c[0]) in couplings_dict[(k+1)%2]:
                    overlaps[frozenset(c)] = (couplings_dict[k][c]+couplings_dict[(k+1)%2][(c[1],c[0])])/2
        overlaps_ordered = OrderedDict({tuple(k): v for k, v in sorted(overlaps.items(), key=lambda item: item)}) # order by mean score
    return exclus+[overlaps_ordered]


def get_normalized_ordered_dict(od):
    fact=1/sum(od.values())
    new_od = OrderedDict()
    for key in od:
        new_od[key] = od[key]*fact
    return new_od
