from comutils.util import *

from urllib.request import urlopen

from Bio import PDB


d_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def get_3to1(three):
    if three in d_3to1:
        return d_3to1[three]
    else:
        raise Exception('Not an amino acid')

def fetch_pdb_file(pdb_id, outputfname):
    try:
        url = "https://files.rcsb.org/download/"+pdb_id+".pdb"
        pdbfile = urlopen(url)
        with open(outputfname+".pdb",'wb') as output:
            output.write(pdbfile.read())
        return outputfname+".pdb"
    except Exception as e:
        url = "https://files.rcsb.org/download/"+pdb_id+".cif"
        ciffile = urlopen(url)
        with open(outputfname+".cif", 'wb') as output:
            output.write(ciffile.read())
        return outputfname+".cif"

def get_pdb_chain(pdb_file, pdb_id, chain_id):
    pdbfile = str(pdb_file)
    if pdbfile.endswith(".pdb"):
        structure = PDB.PDBParser().get_structure(pdb_id, pdbfile)
    elif pdbfile.endswith(".cif"):
        structure = Bio.PDB.MMCIFParser().get_structure(pdb_id, pdbfile)
    else:
        raise Exception("Unknown PDB file format")
    model = structure[0]
    pdb_chain = model[chain_id]
    return pdb_chain

def is_acceptable_residue(residue):
    return residue.get_full_id()[3][0]==' '

def get_res_id(residue):
    return residue.get_full_id()[3][1]

def get_res_letter(residue):
    return get_3to1(residue.get_resname())


def get_seq_pos_to_pdb_chain_pos(sequence, pdb_chain):
    pdb_sequence_dict = {}
    for residue in pdb_chain:
        if is_acceptable_residue(residue):
            pdb_sequence_dict[get_res_id(residue)] = get_res_letter(residue)
    pdb_sequence = "".join(pdb_sequence_dict.values())
    seq_to_pdb_seq = get_pos_first_seq_to_second_seq(sequence, pdb_sequence)
    seq_to_pdb_chain = []
    for seq_index in range(len(sequence)):
        pdb_seq_index = seq_to_pdb_seq[seq_index]
        if pdb_seq_index is None:
            pdb_chain_index = None
        else:
            pdb_chain_index = list(pdb_sequence_dict.keys())[pdb_seq_index]
            assert(get_res_letter(pdb_chain[pdb_chain_index])==sequence[seq_index])
        seq_to_pdb_chain.append(pdb_chain_index)
    return seq_to_pdb_chain


def get_mrf_pos_to_pdb_chain_pos(mrf_pos_to_seq_pos, sequence, pdb_chain):
    mrf_pos_to_pdb_chain_pos = []
    seq_pos_to_pdb_chain_pos = get_seq_pos_to_pdb_chain_pos(sequence, pdb_chain)
    for mrf_pos in range(len(mrf_pos_to_seq_pos)):
        seq_pos = mrf_pos_to_seq_pos[mrf_pos]
        if seq_pos is None:
            pdb_pos = None
        else:
            pdb_pos = seq_pos_to_pdb_chain_pos[seq_pos]
        mrf_pos_to_pdb_chain_pos.append(pdb_pos)
    return mrf_pos_to_pdb_chain_pos


def aa_distance(pos1, pos2, pdb_chain):
    if (pos1 not in pdb_chain) or (pos2 not in pdb_chain):
        raise Exception("Position not in PDB chain")
    r1 = pdb_chain[pos1]
    r2 = pdb_chain[pos2]
    diff_vector = r1['CA'].coord - r2['CA'].coord
    return np.sqrt(np.sum(diff_vector*diff_vector))


def is_pdb_pair_contact(i_pdb, j_pdb, pdb_chain, contact_threshold=8):
    if (i_pdb is None) or (j_pdb is None):
        raise Exception("Position not in PDB chain")
    else:
        return aa_distance(i_pdb, j_pdb, pdb_chain)<=contact_threshold


def is_coupling_contact(i, j, pdb_chain, mrf_pos_to_pdb_chain_pos, contact_threshold=8):
    i_pdb = mrf_pos_to_pdb_chain_pos[i]
    j_pdb = mrf_pos_to_pdb_chain_pos[j]
    return is_pdb_pair_contact(i_pdb, j_pdb, pdb_chain, contact_threshold=contact_threshold)
