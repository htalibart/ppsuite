import pymol
import argparse
import sys

from vizcontacts.contacts_management import *
from comutils import files_management as fm
from comfeature.comfeature import *

def launch_pymol(pdb_id, pdbfile=None):
    pymol.finish_launching(['pymol'])
    if pdbfile is not None:
        pymol.cmd.load(pdbfile)
    else:
        pymol.cmd.fetch(pdb_id, async=0, type="pdb")
    pymol.cmd.show('cartoon')
    pymol.cmd.hide('lines')
    pymol.cmd.hide('nonbonded')
    pymol.cmd.set('cartoon_color', 'grey')
    pymol.cmd.set('cartoon_transparency', 0.2)


def show_coupling(coupling, strength, color, chain_id='A'):
    pos1 = coupling[0]
    pos2 = coupling[1]
    pymol.cmd.select("coupling", "resi "+str(pos1)+" and name CA and chain "+chain_id+" + resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.bond("resi "+str(pos1)+" and name CA and chain "+chain_id, "resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.set_bond("stick_color", color, "coupling")
    pymol.cmd.set_bond("stick_radius", strength, "coupling")
    pymol.cmd.show('sticks', "coupling")
    pymol.cmd.label("coupling", 'resi')


def show_n_couplings(nb_couplings, pdb_seq_couplings_dict, pdb_file, pdb_id, chain_id='A'):
    pdb_chain = fm.get_pdb_chain(pdb_id, pdb_file, chain_id)
    colors = {True : 'blue', False : 'red'}
    for i, (c_set, score) in enumerate(pdb_seq_couplings_dict.items()):
        if i<nb_couplings:
            c = tuple(c_set)
            strength = score
            show_coupling(c, strength, colors[is_true_contact(c, pdb_chain)], chain_id)


def show_predicted_contacts_with_pymol(feature_folder, pdb_id, chain_id='A', nb_couplings=50, pdb_file=None, **kwargs):
    comfeature = ComFeature.from_folder(feature_folder)
    if pdb_file is None:
        pdb_file = fm.get_pdb_file_from_folder(feature_folder)
    pdb_chain = fm.get_pdb_chain(pdb_id, pdb_file, chain_id)
    couplings_dict = get_contact_scores_for_sequence(comfeature)
    pdb_couplings_dict = translate_dict_to_pdb_pos(couplings_dict, pdb_chain, comfeature.sequence)
    launch_pymol(pdb_id, pdb_file)
    show_n_couplings(nb_couplings, pdb_couplings_dict, pdb_file, pdb_id, chain_id=chain_id)


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--feature_folder', help="Feature folder", type=pathlib.Path)
    parser.add_argument('-id', '--pdb_id', help="PDB id")
    parser.add_argument('--pdb_file', help="PDB file")
    args = vars(parser.parse_args(args))

    show_predicted_contacts_with_pymol(**args)

if __name__=="__main__":
    main()

