#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import csv
from makepotts.potts_model import *

def adabmdca_to_ccmpredpy(adabmdca_file, output_filename=None):
    """inputs Potts model weights in adabmDCA format and outputs a Potts model object in CCMpredPy format. Prints result in @output_file if not None"""
    with open(str(adabmdca_file), 'r') as adaf:
        csv_reader = csv.reader(adaf, delimiter=' ')
        J_lines = []
        h_lines = []
        for line in csv_reader:
            if line[0]=='J':
                J_lines.append([int(ind) for ind in line[1:5]]+[float(line[5])])
            else:
                h_lines.append([int(ind) for ind in line[1:3]]+[float(line[3])])
        ncol = len(h_lines)//21
        v = np.zeros((ncol, 21))
        for h_line in h_lines:
            v[h_line[0]][h_line[1]] = h_line[2]
        w = np.zeros((ncol, ncol, 21, 21))
        for J_line in J_lines:
            w[J_line[0],J_line[1],J_line[2],J_line[3]]=J_line[4]
        potts_model = Potts_Model.from_parameters(v, w)


        if output_filename is not None:
            potts_model.to_msgpack(output_filename)

    return potts_model



def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'potts_adabmdca', help='Potts model file in adabmDCA format',
        type=pathlib.Path
    )
    parser.add_argument(
        'potts_ccmpredpy', help='Potts model file in CCMpredPy (and PPalign) format',
        type=pathlib.Path
    )
    args = vars(parser.parse_args(args))

    
    adabmdca_to_ccmpredpy(args['potts_adabmdca'], output_filename=args["potts_ccmpredpy"])


if __name__ == '__main__':
    main()
