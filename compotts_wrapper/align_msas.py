import argparse
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

import files_management as fm

def get_msas_aligned(aligned_positions, train_msa_files, output_msa_file):
    print("merging "+train_msa_files[0]+" and "+train_msa_files[1])
    aligns = [SeqIO.parse(f, "fasta") for f in train_msa_files]
    records = []
    for k in range(2):
        align = aligns[k]
        for record in align:
            new_record = record
            new_seq=""
            for pos in list(aligned_positions.values())[k]:
               new_seq+=record[pos]
            new_record.seq=Seq(new_seq)
            records.append(new_record)
    new_alignment = MultipleSeqAlignment(records)
    with open(output_msa_file, 'w') as f:
        AlignIO.write(new_alignment, f, "fasta")
    print("output can be found at "+output_msa_file)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a1', '--msa_1', help="MSA 1")
    parser.add_argument('-a2', '--msa_2', help="MSA 2")
    parser.add_argument('-o', '--output_msa', help="output MSA")
    parser.add_argument('-a', '--aln_res_file', help="output ComPotts aln file")
    args = vars(parser.parse_args())

    aligned_positions = fm.get_aligned_positions_dict_from_compotts_output_file(args["aln_res_file"])
    get_msas_aligned(aligned_positions, [args["msa_1"], args["msa_2"]], args["output_msa"])
