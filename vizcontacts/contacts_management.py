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
from vizcontacts.pdb_utils import *

def get_contact_scores_for_aln_train_via_ccmpredpy(potts_object):
    """ returns an ordered dictionary of contact scores for positions in the train MSA by calling CCMpredPy (deprecated)"""
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


def get_contact_scores_for_potts_model(pm):
    """ returns an ordered dictionary of contact scores for positions in the train MSA """
    pm.w = pm.w[:,:,:20,:20]
    w_norms = pm.get_w_norms()
    top_ind1, top_ind2 = top_couplings.get_top_pairs(w_norms, reverse_order=False)
    contact_scores = OrderedDict()
    for ind1, ind2 in zip(top_ind1, top_ind2):
        contact_scores[(ind1, ind2)] = w_norms[ind1, ind2]
    return contact_scores




def get_contact_scores_with_pdb_indexes(potts_object, pdb_chain):
    mrf_pos_to_pdb_pos = get_mrf_pos_to_pdb_chain_pos(potts_object.mrf_pos_to_seq_pos, potts_object.sequence, pdb_chain)
    pm_scores = get_contact_scores_for_potts_model(potts_object.potts_model)
    pdb_contact_scores = OrderedDict()
    for c in pm_scores:
        c_pdb = (mrf_pos_to_pdb_pos[c[0]], mrf_pos_to_pdb_pos[c[1]])
        if (c_pdb[0] is not None) and (c_pdb[1] is not None):
            pdb_contact_scores[c_pdb] = pm_scores[c]
    return pdb_contact_scores



#def get_contact_scores_for_sequence(potts_object):
#    """ returns an ordered dictionary of contact scores for positions in the original sequence """
#    aln_scores = get_contact_scores_for_potts_modelpotts_object.potts_model)
#    seq_contact_scores = OrderedDict()
#    for c in aln_scores:
#        c_seq = (potts_object.mrf_pos_to_seq_pos[c[0]], potts_object.mrf_pos_to_seq_pos[c[1]])
#        seq_contact_scores[c_seq] = aln_scores[c]
#    return seq_contact_scores



def get_colored_true_false_dicts(pdb_couplings_dict, pdb_chain, real_sequence, colors={True:'blue', False:'red'}, contact_threshold=8):
    """ tf_d['blue'][c] is the strength of coupling c predicted as true contact and tf_d['red'][c] is the strength of coupling c not predicted as contact """
    tf_d = {colors[val]:OrderedDict() for val in colors}
    for pdb_c in pdb_couplings_dict:
        if not None in pdb_c:
            tf_d[colors[is_true_contact(pdb_c, pdb_chain, contact_threshold=contact_threshold)]][pdb_c] = pdb_couplings_dict[pdb_c]
    return tf_d

def remove_couplings_too_close(couplings_dict, coupling_sep_min):
    """ returns dictionary where couplings separated by less than @coupling_sep_min are removed """
    ok_dict = OrderedDict()
    for c in couplings_dict:
        if None not in c:
            if abs(c[0]-c[1])>coupling_sep_min:
                ok_dict[c] = couplings_dict[c]
    return ok_dict

def get_smaller_dict(couplings_dict, nb_couplings):
    """ returns couplings dictionary with only the strongest @nb_couplings couplings """
    new_dict = OrderedDict()
    for c in couplings_dict:
        if len(new_dict)<nb_couplings:
            new_dict[c] = couplings_dict[c]
    return new_dict


def get_elbow_index(couplings_dict, plot_elbow=False):
    """ find index at which coupling strengths starts to dramatically decrease """
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
    """ returns index at which coupling strength starts to be lower than @score_cutoff in ordered @couplings_dict """
    y = list(couplings_dict.values())
    if y[0]<score_cutoff:
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
