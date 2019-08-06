import numpy as np

from Bio import pairwise2

from basic_modules.global_variables import ALPHABET
q = len(ALPHABET)


def code(c):
        """
            gives number code in [0,21] for letter c
        """
        return ALPHABET.find(c.upper())


def code_whole_seq(sequence):
    return [code(c) for c in sequence]


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


def euclidean_norm(vector):
    return np.linalg.norm(vector)


def scalar_product(v1, v2):
    return np.dot(v1.flatten(), v2.flatten())


def sign_ind(x):
    if (x>=0):
        return 1
    else:
        return -1


def seq_identity(seq1, seq2):
    """ calcule l'identité de séquence entre @seq1 et @seq2 """
    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    nb_matches = alignment[2]
    length_alignment = len(alignment[0])
    #score = (nb_matches/max(len(seq1),len(seq2)))
    score = nb_matches/length_alignment
    return score



def get_trimmed_sequence_for_msa(msa_file, seq):
    """ retourne la séquence @seq taillée pour qu'elle rentre dans le MSA @msa_file : on enlève toutes les positions insérées """
    temp_folder = create_folder("temp_trimming_folder_"+seq)
    ungapped_seq = re.sub('-','',seq) # la séquence sert de nom pour les fichiers, pour éviter les conflits
    tempseq_file = temp_folder+ungapped_seq+".fasta"
    seq_name = create_seq_fasta(seq, tempseq_file)
    combined_aln_file = temp_folder+ungapped_seq+"_aln.fasta"
    os.system("muscle -profile -in1 "+msa_file+" -in2 "+tempseq_file+" -out "+combined_aln_file+" -gapopen -0.1")
    comb_msa = AlignIO.read(combined_aln_file, "fasta")
    new = comb_msa[len(comb_msa)-1].seq
    trimmed=""
    n = len(comb_msa[0].seq)
    for i in range(n):
        if not is_gap_column(i, comb_msa[0:len(comb_msa)-1]): # si pas insertion
            trimmed+=new[i]
    os.system("rm -r "+temp_folder)
    return trimmed



def remove_bad_sequences(input_file, output_file, gap_threshold):
    alignment = AlignIO.read(open(input_file), "fasta")
    acceptable_records = []
    for record in alignment:
        nb_gaps = str(record.seq).count('-')
        if nb_gaps/len(str(record.seq))<1-gap_threshold:
            acceptable_records.append(record)
    with open(output_file, 'w') as f:
        SeqIO.write(acceptable_records, f, "fasta")
