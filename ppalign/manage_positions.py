import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

import comutils.files_management as fm
from comutils.util import *


def aln_dict_to_tuples_list(aln_dict):
    """ input: dictonary of aligned positions in the form {"pos_ref":[,], "pos_2":[,]}, output: tuples list: [(,),...] """
    tuples_list = []
    for pos_aln in range(len(aln_dict["pos_ref"])):
        tuples_list.append((aln_dict["pos_ref"][pos_aln], aln_dict["pos_2"][pos_aln]))
    return tuples_list


def tuples_list_to_aln_dict(tuples_list):
    """ input: tuples list: [(,),...], output: dictonary of aligned positions in the form {"pos_ref":[,], "pos_2":[,]} """
    return {"pos_ref":[pair[0] for pair in tuples_list], "pos_2":[pair[1] for pair in tuples_list]}


def get_alignment_with_gaps(aligned_positions):
    """ input: dictionary of aligned positions, output: dictionary of aligned positions including the gap symbol """
    tuples_list = aln_dict_to_tuples_list(aligned_positions)
    tuples_list_with_gaps = []
    previous_pair = (-1,-1)
    for pair in tuples_list:
        for pos_in_ref in list(range(previous_pair[0]+1,pair[0])):
            tuples_list_with_gaps.append((pos_in_ref,'-'))
        for pos_in_2 in list(range(previous_pair[1]+1,pair[1])):
            tuples_list_with_gaps.append(('-',pos_in_2))
        tuples_list_with_gaps.append(pair)
        previous_pair=pair
    return tuples_list_to_aln_dict(tuples_list_with_gaps)


def get_seq_positions_from_aln_dict(aligned_positions, objects):
    """ input : dict of lists of positions aligned by solver, output : dict of lists of positions in the original sequence """
    real_seq_positions = {}
    c_names = ['pos_ref', 'pos_2']
    for k in range(2):
        real_seq_positions[c_names[k]] = objects[k].get_seq_positions(aligned_positions[c_names[k]])
    real_seq_positions_without_none = {'pos_ref':[], 'pos_2':[]}
    for ind in range(len(real_seq_positions["pos_ref"])):
        if (real_seq_positions["pos_ref"][ind] is not None) and (real_seq_positions["pos_2"][ind] is not None):
            for name in c_names:
                real_seq_positions_without_none[name].append(real_seq_positions[name][ind])
    return real_seq_positions_without_none


def get_aln_positions_from_aln_dict(aligned_positions, objects):
    """ input : dict of lists of positions aligned by solver, output : dict of lists of positions in the original multiple sequence alignment """
    real_aln_positions = {}
    c_names = ['pos_ref', 'pos_2']
    for k in range(2):
        real_aln_positions[c_names[k]] = objects[k].get_aln_positions(aligned_positions[c_names[k]])
    real_aln_positions_without_none = {'pos_ref':[], 'pos_2':[]}
    for ind in range(len(real_aln_positions["pos_ref"])):
        if (real_aln_positions["pos_ref"][ind] is not None) and (real_aln_positions["pos_2"][ind] is not None):
            for name in c_names:
                real_aln_positions_without_none[name].append(real_aln_positions[name][ind])
    return real_aln_positions_without_none



def aligned_positions_to_aligned_sequences(seq_positions, sequences):
    """ translates dict of aligned positions to dict of aligned sequences """
    c_names = ["pos_ref", "pos_2"]
    seqs_aligned = ["",""]
    for k in range(2):
        ck = c_names[k]
        for pos in seq_positions[ck]:
            if (pos=='-'):
                car='-'
            else:
                car=sequences[k][pos]
            seqs_aligned[k]+=car
   # if (len(sequences[0])-seq_positions['pos_ref'][-1]+1>0):
   #     for pos_in_ref in range(seq_positions['pos_ref'][-1]+1,len(sequences[0])):
   #         seqs_aligned[0]+=sequences[0][pos_in_ref]
   #         seqs_aligned[1]+='-'
   # if (len(sequences[1])-seq_positions['pos_2'][-1]+1>0):
   #     for pos_in_2 in range(seq_positions['pos_2'][-1]+1, len(sequences[1])):
   #         seqs_aligned[1]+=sequences[1][pos_in_2]
   #         seqs_aligned[0]+='-'
    return seqs_aligned


def get_seq_positions(aligned_positions, objects):
    """ aligned positions in the models to aligned positions in the sequences """
    seq_positions = {}
    seq_positions['pos_ref'] = objects[0].get_seq_positions(aligned_positions['pos_ref'])
    seq_positions['pos_2'] = objects[1].get_seq_positions(aligned_positions['pos_2'])
    return seq_positions


def remove_None_positions(positions_dict):
    """ remove aligned positions if one of them is None (e.g. if position is in the MSA but not in the sequence) """
    new_positions_dict = {'pos_ref':[], 'pos_2':[]}
    for pos_aln in range(len(positions_dict['pos_ref'])):
        if (positions_dict['pos_ref'][pos_aln] is not None) and (positions_dict['pos_2'][pos_aln]) is not None:
            new_positions_dict['pos_ref'].append(positions_dict['pos_ref'][pos_aln])
            new_positions_dict['pos_2'].append(positions_dict['pos_2'][pos_aln])
    return new_positions_dict


def get_seqs_aligned_using_aln(aligned_positions, objects):
    """ aligned positions in the model to aligned sequences """
    seq_positions_from_objects = get_seq_positions(aligned_positions, objects)
    seq_positions_without_None = remove_None_positions(seq_positions_from_objects)
    seq_positions = get_alignment_with_gaps(seq_positions_without_None)
    return aligned_positions_to_aligned_sequences(seq_positions, [obj.sequence for obj in objects]) 

            
def get_seqs_aligned_in_fasta_file_using_aln(aligned_positions, objects, output_file):
    """ (positions aligned by solver + objects) -> sequences aligned -> in output_file """
    seqs_aligned = get_seqs_aligned_using_aln(aligned_positions, objects)
    seq_records = [SeqRecord(Seq(s), id=o.get_name(), description='') for s,o in zip(seqs_aligned, objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))


def get_seqs_aligned(aligned_positions_with_gaps, objects):
    """ aligned positions in the model with gaps to aligned sequences """
    seq_positions_from_objects = get_seq_positions(aligned_positions_with_gaps, objects)
    seq_positions_without_None = remove_None_positions(seq_positions_from_objects)
    return aligned_positions_to_aligned_sequences(seq_positions_without_None, [obj.sequence for obj in objects]) 


def get_seqs_aligned_in_fasta_file(aligned_positions_with_gaps, objects, output_file):
    """ (positions aligned by solver with gaps + objects) -> sequences aligned -> in output_file """
    seqs_aligned = get_seqs_aligned(aligned_positions_with_gaps, objects)
    seq_records = [SeqRecord(Seq(s), id=o.get_name(), description='') for s,o in zip(seqs_aligned, objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))


def aln_csv_to_aln_fasta(aln_file, objects, output_file):
    """ @aln_file: alignment file provided by PPalign + @objects (Potts objects aligned) -> output fasta file of sequences aligned """
    get_seqs_aligned_in_fasta_file(fm.get_aligned_positions_dict_from_ppalign_output_file(aln_file), objects, output_file)


def aln_sequences_csv_to_aln_fasta(aln_sequences_file, objects, output_file):
    """ @aln_file: csv file of positions in sequences aligned + @objects (Potts objects aligned) -> output fasta file of sequences aligned """
    aligned_positions = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_sequences_file)
    alignment_with_gaps = get_alignment_with_gaps(aligned_positions)
    aligned_sequences = aligned_positions_to_aligned_sequences(alignment_with_gaps, [obj.sequence for obj in objects])
    seq_records = [SeqRecord(Seq(s), id=o.get_name(), description='') for s,o in zip(aligned_sequences, objects)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))

def aln_sequences_csv_to_aln_fasta_using_sequences_only(aln_sequences_file, sequence_files, output_file):
    """ @aln_file: csv file of positions in sequences aligned + @sequences (fasta file for the two sequences aligned) -> output fasta file of sequences aligned """
    aligned_positions = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_sequences_file)
    aligned_positions_without_None = remove_None_positions(aligned_positions)
    alignment_with_gaps = get_alignment_with_gaps(aligned_positions_without_None)
    records = [list(SeqIO.parse(seq_file, 'fasta'))[0] for seq_file in sequence_files]
    aligned_sequences = aligned_positions_to_aligned_sequences(alignment_with_gaps, [str(rec.seq) for rec in records])
    seq_records = [SeqRecord(Seq(s), id=rec.id, description='') for s,rec in zip(aligned_sequences, records)]
    with open(str(output_file), 'w') as f:
        SeqIO.write(seq_records, f, "fasta")
    print("output can be found at "+str(output_file))



def get_pos_aligned_at_pos(aligned_positions_dict, pos):
    """ returns the position in object 2 aligned at position @pos in object 1"""
    d = aligned_positions_dict
    if pos in d['pos_ref']:
        return d['pos_2'][d['pos_ref'].index(pos)]
    else:
        print(str(pos)+" is not aligned")
        return None


def get_initial_positions(aligned_positions, mrf_pos_to_initial_pos_dict):
    """ inputs @aligned_positions: dict of aligned positions, @mrf_pos_to_initial_pos_dict: list of positions in the original sequence for each position in the Potts model
        outputs dict of aligned positions in the original sequences """
    initial_positions = {}
    for key in aligned_positions:
        initial_positions[key]=[]
        for pos in aligned_positions[key]:
            if pos=='-':
                init_pos='-'
            else:
                init_pos = mrf_pos_to_initial_pos_dict[key][pos]
            initial_positions[key].append(init_pos)
    return initial_positions


def get_mrf_pos_to_seq_pos(original_first_seq, seq, mrf_pos_to_aln_pos):
    """ @original_first_seq: typically first sequence in the seed MSA
        @seq: sequence
        outputs list of positions in the sequence for each position in the seed MSA """
    seq_aln_pos = get_pos_first_seq_to_second_seq(original_first_seq, seq)
    mrf_pos_to_seq_pos = [seq_aln_pos[pos] for pos in mrf_pos_to_aln_pos]
    return mrf_pos_to_seq_pos


def get_original_msas_aligned_from_aligned_positions(aligned_positions_dict, objects, output_msa_file):
    """ inputs
        @aligned_positions_dict: dict of aligned positions outputted by PPalign
        @objects: Potts objects that were aligned
        outputs an alignment of the two original MSAs in @output_msa_file """
    msa_aligned_positions_dict = get_aln_positions_from_aln_dict(aligned_positions_dict, objects)
    aligns = [SeqIO.parse(str(obj.seed_aln), "fasta") for obj in objects]
    records = []
    for k, c_name in zip(range(2), ['pos_ref', 'pos_2']):
        align = aligns[k]
        for record in align:
            new_record = record
            new_seq=""
            for pos in msa_aligned_positions_dict[c_name]:
                new_seq+=record[pos]
            new_record.seq=Seq(new_seq)
            records.append(new_record)
    new_alignment = MultipleSeqAlignment(records)
    with open(str(output_msa_file), 'w') as f:
        AlignIO.write(new_alignment, f, "fasta")
    print("output can be found at "+str(output_msa_file))
