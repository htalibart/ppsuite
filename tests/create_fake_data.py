import numpy as np
import random
import os
from basic_modules.global_variables import ALPHABET_WITHOUT_GAP


"""
Permet de créer un faux alignement multiple à partir d'un modèle
Prend en entrée :
- alphabet : un alphabet (dont les lettres serviront à créer l'alignement)
- n : nombre de séquences qu'il y aura dans l'alignement
- proba_noise : probabilité de bruit
- nb_letters_conserved : nombre de lettres différentes par colonne pseudo-conservée (hors bruit)
- template : l'allure de l'alignement multiple. Tableau de chaînes de caractères, chaque position indique le type de colonne :
    - x : colonne totalement aléatoire
    - y : colonne "pseudo-conservée", c'est-à-dire composée d'un ensemble restreint de lettres (de taille nb_letters_conserved)
    - (lettre dans l'alphabet) : colonne conservée, càd avec presque uniquement la lettre en question, et du bruit
    - numéro : colonne couplée avec la colonne correspondant au numéro
    exemple : template=['y','y','6','y','R','9','y','y','11','y','9','y']
- alnfname : nom de fichier de sortie .aln
- nom fichier de sortie .fasta
- nom fichier de sortie .fasta contenant la première séquence du MSA (pour pouvoir appeler les méthodes DCA ensuite)
"""



def one_random_letter(alphabet):
    """ retourne une lettre aléatoire dans l'alphabet donné en entrée """
    return alphabet[random.randint(0,len(alphabet)-1)]

def n_random_letters(alphabet, n):
    """ retourne une liste de n lettres prises aléatoirement dans l'alphabet """
    return random.sample(set(alphabet), n) 

def noise(proba_noise):
    """ retourne vrai si on devrait avoir du bruit, en fonction de la probabilité proba_noise """
    return random.random()<proba_noise


def get_n_majority_letters(column, nb_letters):
    """ retourne les nb_letters lettres majoritaires dans la colonne column """
    counts = {}
    for l in column:
        if l in counts:
            counts[l]+=1
        else:
            counts[l]=1

    majority_letters = []
    for l in sorted(counts, key=counts.get, reverse=True):
        if (len(majority_letters)<nb_letters):
            majority_letters.append(l)
    return majority_letters



def generate_probadict_for_pseudoconserved_column(conserved_letters, proba_noise, alphabet):
    """ génère un dictionnaire de probabilités de type {lettre : proba, lettre : proba, ...} """
    probadict = {}
    proba_not_noise = 1-proba_noise
    proba_each_conserved_letter = proba_not_noise/len(conserved_letters)
    for l in alphabet:
        probadict[l] = proba_noise
    for l in conserved_letters:
        probadict[l] = proba_each_conserved_letter
    return probadict



def generate_column_from_probadict(n, probadict, proba_noise, alphabet):
    """ génère une colonne à partir d'un dictionnaire de probabilités {lettre : proba, lettre : proba...} """
    c = []
    items = probadict.items()
    letters = [item[0] for item in items]
    weights = [item[1] for item in items]
    for k in range(n):
        if noise(proba_noise):
            l = one_random_letter(alphabet)
        else:
            l = np.random.choice(letters, weights)[0]
        c.append(l)
    return c


def get_probadict_from_column(column):
    probadict = {}
    for l in column:
        if l in probadict:
            probadict[l]+=1
        else:
            probadict[l]=1
    for l in probadict:
        probadict[l]=probadict[l]/len(column)
    return probadict

def conserved_column(conserved_letter, n, proba_noise, alphabet):
    """ retourne une colonne conservée (une lettre conserved_letter) de taille n sur un alphabet, avec une probabiltié de bruit proba_noise """
    probadict = generate_probadict_for_pseudoconserved_column([conserved_letter], proba_noise, alphabet)
    return generate_column_from_probadict(n, probadict, proba_noise, alphabet)


def random_column(n, alphabet, proba_noise):
    """ prend un alphabet et une probabilité de bruit prova_noise et retourne une colonne totalement aléatoire de taille n """
    return pseudoconserved_column(alphabet, n, proba_noise, alphabet)


def generate_probadict_for_one_coupling(letters1, letters2):
    """ prend l'alphabet de la première colonne et génère la liste des probas de passer d'une lettre à l'autre, sous forme de dictionnaire {lettre : {{lettre : proba}, {lettre : proba},...}, lettre :...}"""
    probadict = {}
    for l1 in letters1:
        probadict[l1] = {}
        for l2 in letters2:
            probadict[l1][l2] = 1/len(letters2)
    return probadict


def get_correlated_column(c1, n, alphabet, proba_noise, probadict):
    """ prend une colonne c1 d'alphabet letters 1 et retourne une colonne couplée selon le dictionnaire proba dict"""
    c2 = []
    letters1 = list(probadict.keys())
    for l1 in c1:
        if l1 in letters1:
            items = probadict[l1].items()
            letters2 = [item[0] for item in items]
            weights = [item[1] for item in items]
            l2 = random.choices(letters2, weights)[0]
        else:
            l2 = one_random_letter(alphabet)
        c2.append(l2)
    return c2




def print_result_in_file(msa, n, fname, reverse):
    """ prend le MSA créé et l'écrit dans le fichier de taille fname, en inversant l'ordre des colonnes si reverse """
    p = len(msa[0])
    with open(fname, 'w') as f:
        for i in range(p):
            if not reverse:
                for column in msa:
                    f.write(column[i])
                f.write("\n")
            else:
                for column in reversed(msa):
                    f.write(column[i])
                f.write("\n")



def get_aligned_column(ref, n, alphabet, proba_noise, nb_letters, msa, aln_dict, no_coupling):
    """ retourne la colonne alignée à la colonne de référence (dont les infos sont dans @ref) dans l'alignement de référence, en fonction de l'alphabet, des probabilités de mutation proba_noise, de l'alignement en cours, et du dictionnaire aln_dict donnant les correspondances "colonne de l'alignement en cours" : "colonne dans l'alignement de référence
    ref['type'] = type de la colonne
    ref['col'] = la colonne
    ref['ind_col1'] = si coupling, indice de la première colonne du coupling (dans le premier template)
    ref['probadict_coupling'] = si coupling, dictionnaire de probabilités de passage d'une lettre à l'autre : {lettre : {lettre : proba}, {lettre : proba}, ...}
    ref['probadict_column'] = dictionnaire de probabilités de la colonne {lettre : proba, lettre : proba, ...}
 """
    if ref['type']=="x": # si la colonne de référence est de type x (random), on retourne une colonne random
        return random_column(n, alphabet, proba_noise)
    elif ref['type']=="y": # si la colonnne est de type y (pseudo-conservée), on regarde les nb_letters majoritaires de la colonne de référence et on crée une colonne pseudo-conservée à partir d'elles 
        probadict = ref['probadict_column']
        return generate_column_from_probadict(n, probadict, proba_noise, alphabet)
    elif ref['type']=="C": # si la colonne est de type C (conservée), on regarde la lettre la plus conservée dans la référence et on crée une colonne conservée à partir d'elle
        letter = get_n_majority_letters(ref['col'], 1)
        return conserved_column(letter, n, proba_noise, alphabet)
    elif ref['type']=="c2": # si la colonne est la deuxième colonne d'un coupling
        if (no_coupling):    
            probadict = get_probadict_from_column(ref['col'])
            return generate_column_from_probadict(n, probadict, proba_noise, alphabet)
        else:
            ind_col1 = aln_dict[ref['ind_col1']] # on regarde quel est l'indice de la première colonne dans notre alignement
            col1 = msa[ind_col1]
            probadict = ref['probadict_coupling']
            return get_correlated_column(col1, n, alphabet, proba_noise, probadict)


def create_MSA(template, n, alphabet, proba_noise, nb_letters_conserved, reference_dict):
    p = len(template)
    msa = [[]]*p
    letters_list = [[]]*p
    next_list=[]
    previous_list=[]
    new_ref_dict = {}
    aln_dict = {}

    for i in range(p):
        new_ref_dict[i] = {}
        car = template[i]
        if (car=="x"): # colonne totalement random
            msa[i] = random_column(n, alphabet, proba_noise)
            new_ref_dict[i]['type'] = "x"


        elif (car=="y"): # colonne pseudo-conservée
            letters_list[i] = n_random_letters(alphabet, nb_letters_conserved)
            probadict = generate_probadict_for_pseudoconserved_column(letters_list[i], proba_noise, alphabet)
            msa[i] = generate_column_from_probadict(n, probadict, proba_noise, alphabet)
            new_ref_dict[i]['type'] = "y"
            new_ref_dict[i]['col'] = msa[i]
            new_ref_dict[i]['probadict_column'] = probadict

        elif (car in alphabet): # colonne conservée
            msa[i] = conserved_column(car, n, proba_noise, alphabet)
            new_ref_dict[i]['type'] = ["C"]


        elif (car[0]=="["): # colonne devant s'aligner au MSA précédemment calculé
            no_coupling = (car[1]=="-")
            col_nb = int(car[1+no_coupling:len(car)-1])
            aln_dict[col_nb] = i # on stocke un dictionnaire {num ancienne colonne : num nouvelle colonne} (pour les couplings)
            msa[i] = get_aligned_column(reference_dict[col_nb], n, alphabet, proba_noise, nb_letters_conserved, msa, aln_dict, no_coupling)
            letters_list[i] = get_n_majority_letters(msa[i], nb_letters_conserved)
            #new_ref_dict[i]['probadict_column'] = reference_dict[col_nb]['probadict_column']


        else: # colonne faisant partie d'un coupling
            letters_list[i]  = n_random_letters(alphabet, nb_letters_conserved)

            if (i in next_list): # colonne couplée avec une colonne déjà créée
                ind = previous_list[next_list.index(i)]
                letters1 = get_n_majority_letters(msa[ind], nb_letters_conserved)
                probadict = generate_probadict_for_one_coupling(letters1, letters_list[i])
                msa[i] = get_correlated_column(msa[ind], n, alphabet, proba_noise, probadict)
                previous_list.remove(ind)
                next_list.remove(i)
                new_ref_dict[i]['type'] = "c2"
                new_ref_dict[i]['col'] = msa[i]
                new_ref_dict[i]['ind_col1'] = ind
                new_ref_dict[i]['probadict_coupling'] = probadict

            else: # colonne pas encore couplée, qui sera couplée ensuite
                probadict = generate_probadict_for_pseudoconserved_column(letters_list[i], proba_noise, alphabet)
                msa[i] = generate_column_from_probadict(n, probadict, proba_noise, alphabet)
                new_ref_dict[i]['type'] = "y"
                new_ref_dict[i]['col'] = msa[i]
                new_ref_dict[i]['probadict_column'] = probadict

            if(int(car)>i): # si la colonne doit être couplée à une autre colonne plus loin
                next_list.append(int(car))
                previous_list.append(i)
    return msa, new_ref_dict


def create_fake_seq_fasta(aln_file, fastaseq_file): # première séquence d'un aln
    with open(aln_file, 'r') as f:
        lines = f.readlines()

    line = lines[0]
    seq = line[0:len(line)-1]

    with open(fastaseq_file, 'w') as of:
        of.write(">0\n")
        of.write(seq+"\n")

def create_seq_fasta(seq, fastaseq_file): # à partir d'une séquence
    with open(fastaseq_file, 'w') as of:
        of.write(">0\n")
        of.write(seq+"\n")


def create_fake_fasta(aln_file, output_fasta):
    with open(aln_file, 'r') as f:
        lines = f.readlines()

    with open(output_fasta, 'w') as of:
        i=0
        for line in lines:
            of.write(">"+str(i)+"\n")
            seq = line[0:len(line)-1]
            of.write(seq+"\n")
            i+=1



def main(templates, alnfnames, fastafnames, alphabet=ALPHABET_WITHOUT_GAP, n=1000, proba_noise=0.01, nb_letters_conserved=4):
    ref = None
    for i in range(len(templates)):
        msa, new_ref = create_MSA(templates[i], n, alphabet, proba_noise , nb_letters_conserved, ref)
        if i==0:
            ref = new_ref
        print_result_in_file(msa, n, alnfnames[i], False)
        create_fake_fasta(alnfnames[i], fastafnames[i])


if __name__=="__main__":
    alphabet = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"];

    n = 1000
    proba_noise = 0.01
    nb_letters_conserved = 4

    ref_template = ["1", "0", "3", "2"]
    template2 = ["[0]", "y", "[1]", "[2]"]
    templates = [ref_template, template2]


    alnfnames = ["fake_aln_1","fake_aln_2"]
    fastafnames = ["fake_fasta_1","fake_fasta_2"]
    mrfnames = ["fake_mrf_1", "fake_mrf_2"]

    main(alphabet, n, proba_noise, nb_letters_conserved, templates, alnfnames, fastafnames, mrfnames)
