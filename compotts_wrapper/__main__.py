import argparse

from call_compotts import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', '--potts_model_1', help="Potts Model 1")
    parser.add_argument('-p2', '--potts_model_2', help="Potts Model 2")
    parser.add_argument('-a', '--aln_res_file', help="Output alignment file")
    parser.add_argument('-i', '--info_res_file', help="Output info file")
    args = vars(parser.parse_args())


    aln_res_file = args["aln_res_file"]
    info_res_file = args["info_res_file"]
    align_two_potts_models_from_files([args["potts_model_1"], args["potts_model_2"]], aln_res_file, info_res_file)
