""" module mainly used for test purposes: create artificial MSA for test """
import numpy as np
import random
import os
from comutils.global_variables import ALPHABET_WITHOUT_GAP


"""
Input :
- alphabet
- n : number of sequences in the output MSA
- proba_noise : probability for noise in each column
- nb_letters_conserved : number of letters in each conserved column (noise aside)
- template: list of characters, on per column in the MSA, indicating the column type, where:
    - 'x' means completely random column (letters are randomly sampled from alphabet)
    - 'y' means conserved column, i.e. there are only @nb_letters_conserved in the column (+ noise)
    - (alphabet letter): conserved letter with only this one letter (+ noise)
    - (number): column coupled with column at position (number)
    example : template=['y','y','6','y','R','9','y','y','11','y','9','y']
- alnfname : output file .aln
- output file .fasta
"""


def one_random_letter(alphabet):
    """ returns a random letter in @alphabet """
    return alphabet[random.randint(0,len(alphabet)-1)]

def n_random_letters(alphabet, n):
    """ returns @n random letters in @alphabet """
    return random.sample(set(alphabet), n) 

def noise(proba_noise):
    """ returns True if this letter should be noise if the probability of noise is @proba_noise """
    return random.random()<proba_noise


def get_n_majority_letters(column, nb_letters):
    """ returns the most frequent @nb_letters letters in @column """
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
    """ outputs a probability dict like {letter:proba, letter:proba, ...} for each letter in the alphabet to build a column with a list @conserved_letters of conserved letters and a probability of noise @proba_noise"""
    probadict = {}
    proba_not_noise = 1-proba_noise
    proba_each_conserved_letter = proba_not_noise/len(conserved_letters)
    proba_each_noise = proba_noise/(len(alphabet)-len(conserved_letters))
    for l in alphabet:
        probadict[l] = proba_each_noise
    for l in conserved_letters:
        probadict[l] = proba_each_conserved_letter
    return probadict



def generate_column_from_probadict(n, probadict, proba_noise, alphabet):
    """ outputs a column from a dictionary providing a probability for each letter in the alphabet """
    c = []
    items = probadict.items()
    letters = [item[0] for item in items]
    weights = [item[1] for item in items]
    for k in range(n):
        if noise(proba_noise):
            l = one_random_letter(alphabet)
        else:
            l = np.random.choice(letters, p=weights)[0]
        c.append(l)
    return c


def get_probadict_from_column(column):
    """ outputs a probability dict {letter:proba, ...} given a column """
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
    """ outputs a column of depth @n with one conserved letter @conserved_letter with noise probability @proba_noise """
    probadict = generate_probadict_for_pseudoconserved_column([conserved_letter], proba_noise, alphabet)
    return generate_column_from_probadict(n, probadict, proba_noise, alphabet)


def random_column(n, alphabet, proba_noise):
    """ outputs a column of depth @n with completely random letters """
    return pseudoconserved_column(alphabet, n, proba_noise, alphabet)


def generate_probadict_for_one_coupling(letters1, letters2):
    """ outputs a dictionary {letter1:{letter2:proba},...} where letter1 is in @letters1 and letter2 is in @letters2 """
    probadict = {}
    for l1 in letters1:
        probadict[l1] = {}
        for l2 in letters2:
            probadict[l1][l2] = 1/len(letters2)
    return probadict


def get_correlated_column(c1, n, alphabet, proba_noise, probadict):
    """ inputs a column @c1 of depth @n on @alphabet and outputs a coupled column according to a probability dictionary @probadict in the form {letter:{letter:proba},...}"""
    c2 = []
    letters1 = list(probadict.keys())
    for l1 in c1:
        if l1 in letters1:
            items = probadict[l1].items()
            letters2 = [item[0] for item in items]
            weights = [item[1] for item in items]
            l2 = np.random.choice(letters2, p=weights)[0]
        else:
            l2 = one_random_letter(alphabet)
        c2.append(l2)
    return c2




def print_result_in_file(msa, n, fname, reverse):
    """ writes created MSA of depth @n in file @fname, reverses order if @reverse """
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
    """ outputs a column in a second MSA that should be aligned to a column in the first MSA described by @ref where:
        - ref['type'] is the column type ('x', 'y', ...)
        - ref['col'] is the actual column
        - ref['ind_col1']: if column is coupled, position of the first coupled column
        - ref['probadict_coupling']: if column is coupled, dict of probabilities {letter:{letter:proba},...}
        - ref['probadict_column']: dict of probabilities for the column {letter:proba, letter:proba, ...}
        @aln_dict is a dictionary providing columns that should be matched, in the form {column in this MSA: column in the first MSA}
        @no_coupling is True if we don't want this column to be coupled even if @ref is
    """
    if ref['type']=="x": # if ref column type is 'x' (random), matched column is random too
        return random_column(n, alphabet, proba_noise)
    elif ref['type']=="y": # if ref column is 'y' (conserved), retrieve the main letters and output column with these letters 
        probadict = ref['probadict_column']
        return generate_column_from_probadict(n, probadict, proba_noise, alphabet)
    elif ref['type']=="C": # if ref letter is 'C' (one conserved letter), output column with the conserved letter
        letter = get_n_majority_letters(ref['col'], 1)
        return conserved_column(letter, n, proba_noise, alphabet)
    elif ref['type']=="c2": # if column is the second column in a coupling
        if (no_coupling):    
            probadict = get_probadict_from_column(ref['col'])
            return generate_column_from_probadict(n, probadict, proba_noise, alphabet)
        else:
            ind_col1 = aln_dict[ref['ind_col1']] # index of the first coupled column in this MSA
            col1 = msa[ind_col1]
            probadict = ref['probadict_coupling']
            return get_correlated_column(col1, n, alphabet, proba_noise, probadict)



def create_MSA(template, n, alphabet, proba_noise, nb_letters_conserved, reference_dict):
    """ creates artificial MSA and outputs a reference dictionary describing it so that the next MSA may be aligned to it 
        @reference_dict gives: {position:{type:, col: probadict_column,...},...} """
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
        if (car=="x"): # completely random conserved
            msa[i] = random_column(n, alphabet, proba_noise)
            new_ref_dict[i]['type'] = "x"


        elif (car=="y"): # conserved column with @nb_letters_conserved letters
            letters_list[i] = n_random_letters(alphabet, nb_letters_conserved)
            probadict = generate_probadict_for_pseudoconserved_column(letters_list[i], proba_noise, alphabet)
            msa[i] = generate_column_from_probadict(n, probadict, proba_noise, alphabet)
            new_ref_dict[i]['type'] = "y"
            new_ref_dict[i]['col'] = msa[i]
            new_ref_dict[i]['probadict_column'] = probadict

        elif (car in alphabet): # conserved column with one letter
            msa[i] = conserved_column(car, n, proba_noise, alphabet)
            new_ref_dict[i]['type'] = ["C"]


        elif (car[0]=="["): # column that will be aligned to a column in a previously built MSA described in @reference_dict
            no_coupling = (car[1]=="-")
            col_nb = int(car[1+no_coupling:len(car)-1])
            aln_dict[col_nb] = i # storing a dictonary {position first column: position second column} for couplings
            msa[i] = get_aligned_column(reference_dict[col_nb], n, alphabet, proba_noise, nb_letters_conserved, msa, aln_dict, no_coupling)
            letters_list[i] = get_n_majority_letters(msa[i], nb_letters_conserved)


        else: # column in a coupling
            letters_list[i]  = n_random_letters(alphabet, nb_letters_conserved)

            if (i in next_list): # column coupled with column already created
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

            else: # first column of a coupling
                probadict = generate_probadict_for_pseudoconserved_column(letters_list[i], proba_noise, alphabet)
                msa[i] = generate_column_from_probadict(n, probadict, proba_noise, alphabet)
                new_ref_dict[i]['type'] = "y"
                new_ref_dict[i]['col'] = msa[i]
                new_ref_dict[i]['probadict_column'] = probadict

            if(int(car)>i): # column will be coupled to a column that will be created later
                next_list.append(int(car))
                previous_list.append(i)

    return msa, new_ref_dict



def create_fake_fasta(aln_file, output_fasta):
    """ converts created aln to fasta """
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
