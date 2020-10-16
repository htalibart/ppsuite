#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import sys
import csv
from Bio import AlignIO
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein


def fasta2indices(alignment_handle, start_at_1=False):
    align = AlignIO.read(alignment_handle, 'fasta')
    tuple_list=[]
    assert(len(align) == 2)
    if start_at_1:
        pos0, pos1 = 0, 0
    else:
        pos0, pos1 = -1, -1
    for pair in zip(align[0].seq, align[1].seq):
        aligned_residues = True
        if pair[0].upper() in ExtendedIUPACProtein.letters:
            pos0 += 1
        else:
            aligned_residues = False
        if pair[1].upper() in ExtendedIUPACProtein.letters:
            pos1 += 1
        else:
            aligned_residues = False
        if aligned_residues:
            tuple_list.append((pos0, pos1))
    return tuple_list


def fasta2csv(alignment_handle, output_handle, start_at_1=False):
    tuple_list = fasta2indices(alignment_handle, start_at_1=start_at_1)
    csv_writer = csv.writer(output_handle)
    csv_writer.writerow(['pos_ref','pos_2'])
    for pair in tuple_list:
        csv_writer.writerow(pair)



def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_alignment', help='Pairwise alignment in fasta format',
        type=argparse.FileType('r')
    )
    parser.add_argument(
        'output_csv', help='Output file in (csv)',
        type=argparse.FileType('w')
    )
    args = vars(parser.parse_args(args))
    fasta2csv(args['input_alignment'], args["output_csv"])


if __name__ == '__main__':
    main()
