import argparse
import sys

from comutils.global_variables import ALPHABET
from vizpotts.vizpotts import *
from comfeature.comfeature import *
from vizcontacts.contacts_management import *

def get_relative_w_norms(comfeature, pdb_chain, contact_distance=8):
    w_norms = comfeature.potts_model.get_w_norms()
    rtpdb = get_real_pos_to_pdb_pos(pdb_chain, comfeature.sequence)
    for i in range(comfeature.potts_model.ncol):
        for j in range(i, comfeature.potts_model.ncol):
            pdb_sequence_coupling = tuple([rtpdb[comfeature.mrf_pos_to_seq_pos[k]] for k in [i,j]])
            if not is_true_contact(pdb_sequence_coupling, pdb_chain, contact_distance):
                w_norms[i][j] = -w_norms[i][j]
                w_norms[j][i] = w_norms[i][j]
    return w_norms


def visualize_mrf_with_contact_map(comfeature, pdb_chain, alphabet=ALPHABET, start_at_1=True, show_figure=True, **kwargs):
    relative_w_norms = get_relative_w_norms(comfeature, pdb_chain)
    visualize_parameters(comfeature.potts_model.v, comfeature.potts_model.get_v_norms(), relative_w_norms, comfeature.get_name(), alphabet=alphabet, start_at_1=start_at_1, show_figure=show_figure)


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--feature_folder', help="Feature folder", type=pathlib.Path)
    parser.add_argument('--pdb_file', help="PDB file", type=pathlib.Path, default=None)
    parser.add_argument('-i', '--pdb_id', help="PDB file")
    parser.add_argument('-cid', '--chain_id', help="PDB chain id", default='A')

    args = vars(parser.parse_args(args))

    comfeature = ComFeature.from_folder(args['feature_folder'])
    if args['pdb_file'] is None:
        name = str(comfeature.folder)+'/'+args['pdb_id']
        args['pdb_file'] = fm.fetch_pdb_file(args['pdb_id'], name)
    pdb_chain = fm.get_pdb_chain(args['pdb_id'], args['pdb_file'], chain_id=args['chain_id'])

    visualize_mrf_with_contact_map(comfeature, pdb_chain, **args)

if __name__=="__main__":
    main()
    parser = argparse.ArgumentParser()
