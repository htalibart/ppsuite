import os
import argparse
import pathlib

import files_management as fm
from multiple_alignment import *


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--seq_folder', help="Folder containing the sequences", type=pathlib.Path)
    parser.add_argument('-a3m', '--a3m_folder', help="Folder containing the corresponding a3m files (must have the same names)", type=pathlib.Path)
    parser.add_argument('-o', '--output_folder', help="Output_folder", type=pathlib.Path)
    # TODO ajouter les autres arguments
    args = vars(parser.parse_args())

    output_folder = args["output_folder"]
    if output_folder is None:
        output_folder = time.strftime("%Y%m%d-%H%M%S")
    fm.create_folder(output_folder)

    seq_folder = args['seq_folder']
    a3m_folder = args['a3m_folder']

    seq_files = [seq_folder+f for f in os.listdir(seq_folder)]
    a3m_files = [a3m_folder+f for f in os.listdir(a3m_folder)]

    multiple_alignment(seq_files, a3m_files, output_folder)
