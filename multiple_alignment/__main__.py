import argparse

import files_management as fm
from multiple_alignment.multiple_alignment import *


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--seq_folder', help="Folder containing the sequences")
    parser.add_argument('-a3m', '--a3m_folder', help="Folder containing the corresponding a3m files (must have the same names)")
    parser.add_argument('-o', '--output_folder', help="Output_folder")
    # TODO ajouter les autres arguments
    args = vars(parser.parse_args())

    output_folder = args["output_folder"]
    if output_folder is None:
        output_folder = time.strftime("%Y%m%d-%H%M%S")
    fm.create_folder(output_folder)

    seq_files = os.listdir(args['seq_folder'])
    a3m_files = os.listdir(args['a3m_folder'])

    multiple_alignment(seq_files, a3m_files, output_folder)
