from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Main
import csv

from comutils import files_management as fm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from makepotts.handle_insertions import *

import pkg_resources
julia_script_insertions_file = pkg_resources.resource_filename('makepotts', 'infer_insertion_penalties.jl')



def infer_insertion_penalties_in_file_using_dcabuild(seed_a3m_file, seed_length, output_file):
    """ calls DCAbuild function to infer insertion penalties from MSA @seed_a3m_file with insertions indicated as lower case letters wrt MSA of length @seed_length, output .tsv file in DCAbuild format @output_file"""
    fm.check_if_file_ok(seed_a3m_file)
    call_dcabuild_infer_ins = Main.include(julia_script_insertions_file)
    call_dcabuild_infer_ins(str(seed_a3m_file), seed_length, str(output_file))
    insertion_penalties = get_insertion_penalties_from_file(output_file)
    # set external insertions to 0
    #for insertion_type in ['open', 'extend']:
    #    insertion_penalties[insertion_type][0] = 0
    #    insertion_penalties[insertion_type].append(0)
    write_insertion_penalties_in_file(insertion_penalties, output_file)

