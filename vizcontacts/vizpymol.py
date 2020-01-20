import argparse
import sys

from vizcontacts.contacts_management import *
from comutils import files_management as fm
from makepotts.potts_object import *

import pymol

def launch_pymol(pdb_id, pdbfile=None):
    #pymol.finish_launching(['pymol'])
    if pdbfile is not None:
        pymol.cmd.load(pdbfile)
    else:
        pymol.cmd.fetch(pdb_id, async_=0, type="pdb")
    pymol.cmd.show('cartoon')
    pymol.cmd.hide('lines')
    pymol.cmd.hide('nonbonded')
    pymol.cmd.set('cartoon_color', 'grey')
    pymol.cmd.set('cartoon_transparency', 0.2)


def show_coupling(pdb_coupling, strength, color, chain_id='A'):
    pos1 = pdb_coupling[0]
    pos2 = pdb_coupling[1]
    pymol.cmd.select("coupling", "resi "+str(pos1)+" and name CA and chain "+chain_id+" + resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.bond("resi "+str(pos1)+" and name CA and chain "+chain_id, "resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.set_bond("stick_color", color, "coupling")
    pymol.cmd.set_bond("stick_radius", strength, "coupling")
    pymol.cmd.show('sticks', "coupling")
    pymol.cmd.label("coupling", 'resi')


def show_n_couplings(nb_couplings, pdb_seq_couplings_dict, pdb_file, pdb_id, chain_id='A', coupling_sep_min=2, thickness=1, colors={True:'green', False:'red'}):
    pdb_chain = fm.get_pdb_chain(pdb_id, pdb_file, chain_id)
    n=0
    for i, (c, score) in enumerate(pdb_seq_couplings_dict.items()):
        if n<nb_couplings:
            if abs(c[0]-c[1])>coupling_sep_min:
                strength = score*thickness
                show_coupling(c, strength, colors[is_true_contact(c, pdb_chain)], chain_id)
                n+=1


def show_predicted_contacts_with_pymol(feature_folders, pdb_id, chain_id='A', pdb_file=None, top=20, coupling_sep_min=3, thickness=1, auto_top=False, wij_cutoff=None, normalize=False, debug_mode=False, **kwargs):

    potts_objects = []
    for feature_folder in feature_folders:
        potts_objects.append(Potts_Object.from_folder(feature_folder))

    if pdb_file is None:
        name = str(potts_objects[0].folder)+'/'+pdb_id
        pdb_file = fm.fetch_pdb_file(pdb_id, name)
    pdb_chain = fm.get_pdb_chain(pdb_id, pdb_file, chain_id)

    pdb_couplings_dicts = []
    tops = []
    for potts_object in potts_objects:
        couplings_dict = get_contact_scores_for_sequence(potts_object)
        pdb_couplings_dict = translate_dict_to_pdb_pos(couplings_dict, pdb_chain, potts_object.sequence)
        if auto_top:
            cutindex = get_elbow_index(pdb_couplings_dict, plot_elbow=debug_mode)
            pdb_couplings_dict = get_smaller_dict(pdb_couplings_dict, cutindex)
            nb_couplings = len(pdb_couplings_dict)
        if wij_cutoff:
            cutindex = get_cutoff_smaller_than(pdb_couplings_dict, wij_cutoff)
            pdb_couplings_dict = get_smaller_dict(pdb_couplings_dict, cutindex)
            nb_couplings = len(pdb_couplings_dict)
        else:
            nb_couplings = top
        pdb_couplings_dict = remove_couplings_too_close(pdb_couplings_dict, coupling_sep_min)
        pdb_couplings_dicts.append(pdb_couplings_dict)
        tops.append(nb_couplings)

    launch_pymol(pdb_id, pdb_file)

    exclus_overlap = get_exclus_overlaps(pdb_couplings_dicts, tops)

    if normalize:
        for k in range(len(exclus_overlap)):
            exclus_overlap[k] = get_normalized_ordered_dict(exclus_overlap[k])

    for d, colors in zip(exclus_overlap, [{True:'green', False:'red'}, {True:'blue', False:'yellow'}, {True:'teal', False:'orange'}]):
        show_n_couplings(len(d), d, pdb_file, pdb_id, chain_id=chain_id, coupling_sep_min=coupling_sep_min, thickness=thickness, colors=colors)


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--feature_folders', help="Feature folder(s)", type=pathlib.Path, nargs='+')
    parser.add_argument('--pdb_file', help="PDB file", type=pathlib.Path, default=None)
    parser.add_argument('-id', '--pdb_id', help="PDB id")
    parser.add_argument('-cid', '--chain_id', help="PDB chain id", default='A')
    parser.add_argument('-sep', '--coupling_sep_min', help="Min. nb residues between members of a coupling", default=3, type=int)
    parser.add_argument('-n', '--top', help="Nb of couplings displayed", type=int, default=20)
    parser.add_argument('--wij_cutoff', help="||wij|| <= wij_cutoff are removed", default=None, type=float) 
    parser.add_argument('--auto_top', help="Nb couplings displayed = elbow of the score curve", default=False, action='store_true')
    parser.add_argument('-t', '--thickness', help="Couplings thickness factor", type=float, default=1)
    parser.add_argument('--normalize', help="Normalize coupling values", default=False, action='store_true')
    parser.add_argument('--debug_mode', help="Debug mode", default=False, action='store_true')
    parser.add_argument('--out_session_file', '-pse', help="PyMOL output session file (must end in .pse)", type=pathlib.Path, default=pathlib.Path('/tmp/tmp_pymol_session_file.pse'))

    args = vars(parser.parse_args(args))

    show_predicted_contacts_with_pymol(**args)
    pymol.cmd.save(str(args["out_session_file"]))
    print("PyMOL session saved at "+str(args["out_session_file"]))

if __name__=="__main__":
    main()

