import pymol

from vizcontacts.contacts_management import *
from comutils import files_management as fm

def launch_pymol(pdbid, pdbfile=None):
    pymol.finish_launching(['pymol'])
    if pdbfile is not None:
        pymol.cmd.load(pdbfile)
    else:
        pymol.cmd.fetch(pdbid, async=0, type="pdb")
    pymol.cmd.show('cartoon')
    pymol.cmd.hide('lines')
    pymol.cmd.hide('nonbonded')
    pymol.cmd.set('cartoon_color', 'grey')
    pymol.cmd.set('cartoon_transparency', 0.2)


def show_coupling(coupling, strength, color, chain_id='A'):
    pos1 = coupling[0]+1
    pos2 = coupling[1]+1
    pymol.cmd.select("coupling", "resi "+str(pos1)+" and name CA and chain "+chain_id+" + resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.bond("resi "+str(pos1)+" and name CA and chain "+chain_id, "resi "+str(pos2)+" and name CA and chain "+chain_id)
    pymol.cmd.set_bond("stick_color", color, "coupling")
    pymol.cmd.set_bond("stick_radius", strength, "coupling")
    pymol.cmd.show('sticks', "coupling")
    pymol.cmd.label("coupling", 'resi')


def show_n_couplings(nb_couplings, couplings_dict, pdbfile, pdbid, chain_id='A'):
    pdb_chain = fm.get_pdb_chain(pdbid, pdbfile, chain_id)
    colors = {True : 'blue', False : 'red'}
    for i, (c_set, score) in enumerate(couplings_dict.items()):
        if i<nb_couplings:
            c = tuple(c_set)
            strength = score
            show_coupling(c, strength, colors[is_true_contact(c, pdb_chain)], chain_id)

