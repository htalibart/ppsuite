import tempfile
import shutil
import numpy as np

from Bio import pairwise2, AlignIO

from comutils.tool_wrapper import *
from comutils import files_management as fm
from comutils.global_variables import ALPHABET
q = len(ALPHABET)


def code(c):
        """
            gives number code in [0,21] for letter c
        """
        return ALPHABET.find(c.upper())


def code_whole_seq(sequence):
    return [code(c) for c in sequence]


def euclidean_norm(vector):
    return np.linalg.norm(vector)


def scalar_product(v1, v2):
    return np.dot(v1.flatten(), v2.flatten())


def sign_ind(x):
    return 2*(x>=0)-1


def seq_identity(seq1, seq2):
    """ calcule l'identité de séquence entre @seq1 et @seq2 """
    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    nb_matches = alignment[2]
    length_alignment = len(alignment[0])
    #score = (nb_matches/max(len(seq1),len(seq2)))
    score = nb_matches/length_alignment
    return score


def is_gap_column(i, msa):
    """ returns True iff the column @i of @msa contains only gaps """
    gap = True
    s=0
    while ( (gap==True) and (s<len(msa)) ):
        seq = msa[s].seq
        if (seq[i]!="-"):
            gap = False
        s+=1
    return gap


def get_trimmed_sequence_for_msa(msa_file, seq):
    """ retourne la séquence @seq taillée pour qu'elle rentre dans le MSA @msa_file : on enlève toutes les positions insérées """
    temp_folder = pathlib.Path(tempfile.mkdtemp())
    tempseq_file = temp_folder/("seq.fasta")
    seq_name = fm.create_seq_fasta(seq, tempseq_file)
    combined_aln_file = temp_folder/"aln.fasta"
    call_muscle_profile(msa_file,tempseq_file,combined_aln_file)
    comb_msa = AlignIO.read(str(combined_aln_file), "fasta")
    new = comb_msa[len(comb_msa)-1].seq
    trimmed=""
    n = len(comb_msa[0].seq)
    for i in range(n):
        if not is_gap_column(i, comb_msa[0:len(comb_msa)-1]): # si pas insertion
            trimmed+=new[i]
    shutil.rmtree(temp_folder)
    return trimmed



def remove_sequences_with_too_many_gaps(input_file, output_file, gap_threshold):
    alignment = AlignIO.read(open(input_file), "fasta")
    acceptable_records = []
    for record in alignment:
        nb_gaps = str(record.seq).count('-')
        if nb_gaps/len(str(record.seq))<1-gap_threshold:
            acceptable_records.append(record)
    with open(output_file, 'w') as f:
        SeqIO.write(acceptable_records, f, "fasta")
