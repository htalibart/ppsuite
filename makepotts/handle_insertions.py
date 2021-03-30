from julia import Main
import csv

from comutils import files_management as fm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


import pkg_resources
julia_script_insertions_file = pkg_resources.resource_filename('makepotts', 'infer_insertion_penalties.jl')

def infer_insertion_penalties_in_file(seed_a3m_file, seed_length, output_file):
    fm.check_if_file_ok(seed_a3m_file)
    call_dcabuild_infer_ins = Main.include(julia_script_insertions_file)
    call_dcabuild_infer_ins(str(seed_a3m_file), seed_length, str(output_file))



def get_insertion_penalties_from_file(insertion_penalties_file):
    fm.check_if_file_ok(insertion_penalties_file)
    insertion_penalties = {'open':[], 'extend':[]}

    # read penalties from file
    with open(str(insertion_penalties_file), 'r') as tsv_file:
        csv_reader = csv.reader(tsv_file, delimiter='\t')
        for row in csv_reader:
            insertion_penalties['open'].append(float(row[0]))
            insertion_penalties['extend'].append(float(row[1]))

    # external gaps:
    for insertion_type in ['open', 'extend']:
        insertion_penalties[insertion_type][0] = 0
        insertion_penalties[insertion_type].append(0)

    return insertion_penalties



def lower_case_trimmed_columns(aln_with_insertions, output_file, list_of_trimmed_columns, fileformat='fasta'):
    records = list(SeqIO.parse(str(aln_with_insertions),fileformat))
    output_records = []
    for record in records:
        upper_index=-1
        new_record = record
        sequence_str = str(record.seq)
        new_sequence_str=''
        for letter in sequence_str:
            if letter.isupper():
                upper_index+=1
                if upper_index not in list_of_trimmed_columns:
                    letter = letter.lower()
            new_sequence_str+=letter
        new_record.seq = Seq(new_sequence_str, IUPAC.protein)
        output_records.append(new_record)
    SeqIO.write(output_records, str(output_file), 'fasta')
