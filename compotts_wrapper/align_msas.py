import argparse
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

def get_msas_aligned(res_aln_file, train_msa_files, output_msa):
    print("merging "+train_msa_files[0]+" and "+train_msa_files[1])
    df = pd.read_csv(res_aln_file)
    aligned_positions = [df['pos_ref'].tolist(), df['pos_2'].tolist()]
    aligns = [SeqIO.parse(f, "fasta") for f in train_msa_files]
    records = []
    for k in range(2):
        align = aligns[k]
        for record in align:
            new_record = record
            new_seq=""
            for pos in aligned_positions[k]:
               new_seq+=record[pos]
            new_record.seq=Seq(new_seq)
            records.append(new_record)
    new_alignment = MultipleSeqAlignment(records)
    with open(output_msa, 'w') as f:
        AlignIO.write(new_alignment, f, "fasta")


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a1', '--msa_1', help="MSA 1")
    parser.add_argument('-a2', '--msa_2', help="MSA 2")
    parser.add_argument('-o', '--output_msa', help="output MSA")
    parser.add_argument('-a', '--aln_res_file', help="output ComPotts aln file")
    args = vars(parser.parse_args())

    get_msas_aligned(args["aln_res_file"], [args["msa_1"], args["msa_2"]], args["output_msa"])
