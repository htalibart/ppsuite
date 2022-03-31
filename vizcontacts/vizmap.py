import argparse
import sys

from comutils.global_variables import ALPHABET
from comutils import files_management as fm
from vizpotts.vizpotts import *
from makepotts.potts_object import *
from vizcontacts.contacts_management import *

def get_relative_w_norms(potts_object, pdb_file, chain_id, contact_distance=8):
    w_norms = potts_object.potts_model.get_w_norms()
    rtpdb = get_real_pos_to_pdb_pos(pdb_file, chain_id, potts_object.sequence)
    for i in range(potts_object.potts_model.ncol):
        for j in range(i, potts_object.potts_model.ncol):
            pdb_sequence_coupling = tuple([rtpdb[potts_object.mrf_pos_to_seq_pos[k]] for k in [i,j]])
            if pdb_sequence_coupling[0] is not None and pdb_sequence_coupling[1] is not None:
                if not is_true_contact(pdb_sequence_coupling, pdb_file, chain_id, contact_distance):
                    w_norms[i][j] = -w_norms[i][j]
                    w_norms[j][i] = w_norms[i][j]
    return w_norms


def visualize_mrf_with_contact_map(potts_object, pdb_file, chain_id, alphabet=ALPHABET, start_at_1=True, show_figure=True, **kwargs):
    relative_w_norms = get_relative_w_norms(potts_object, pdb_file, chain_id)
    visualize_parameters(potts_object.potts_model.v, potts_object.potts_model.get_v_norms(), relative_w_norms, potts_object.get_name(), alphabet=alphabet, start_at_1=start_at_1, show_figure=show_figure)


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--potts_folder', help="Feature folder", type=pathlib.Path, required=True)
    parser.add_argument('--pdb_file', help="PDB file", type=pathlib.Path, default=None)
    parser.add_argument('-i', '--pdb_id', help="PDB file", required=True)
    parser.add_argument('-cid', '--chain_id', help="PDB chain id", default='A')

    args = vars(parser.parse_args(args))

    potts_object = Potts_Object.from_folder(args['potts_folder'])
    if args['pdb_file'] is None:
        name = str(potts_object.folder)+'/'+args['pdb_id']
        args['pdb_file'] = fm.fetch_pdb_file(args['pdb_id'], name)
    fm.check_if_file_ok(args['pdb_file'])
#    pdb_chain = fm.get_pdb_chain(args['pdb_file'], chain_id=args['chain_id'])

    visualize_mrf_with_contact_map(potts_object, **args)

if __name__=="__main__":
    main()
    parser = argparse.ArgumentParser()
