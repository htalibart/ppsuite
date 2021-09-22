#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import csv
from makepotts.potts_model import *

def adabmdca_to_ccmpredpy(adabmdca_file, **kwargs):
    """inputs Potts model weights in adabmDCA format and outputs a Potts model object in CCMpredPy format. Prints result in @output_file if not None"""
    potts_model = Potts_Model.from_adabmdca_file(adabmdca_file, **kwargs)
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

    
    adabmdca_to_ccmpredpy(args['potts_adabmdca'], binary_file=args["potts_ccmpredpy"])


if __name__ == '__main__':
    main()
