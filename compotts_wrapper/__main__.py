import argparse

from call_compotts import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', '--potts_model_1', help="Potts model 1")
    parser.add_argument('-p2', '--potts_model_2', help="Potts model 2")
    parser.add_argument('-s1', '--sequence_file_1', help="Sequence file 1")
    parser.add_argument('-s2', '--sequence_file_2', help="Sequence file 2")
    parser.add_argument('-h1', '--a3m_file_1', help="HH-blits output file 1")
    parser.add_argument('-h2', '--a3m_file_2', help="HH-blits output file 2")
    parser.add_argument('-a', '--aln_res_file', help="Output alignment file")
    parser.add_argument('-i', '--info_res_file', help="Output info file")
    args = vars(parser.parse_args())


    aln_res_file = args["aln_res_file"]
    info_res_file = args["info_res_file"]
    
    if (args["potts_model_1"] is not None) and (args["potts_model_2"] is not None):
        align_two_potts_models_from_files([args["potts_model_1"], args["potts_model_2"]], aln_res_file, info_res_file)
    elif (args["sequence_file_1"] is not None) and (args["sequence_file_2"] is not None) and (args["a3m_file_1"] is not None) and (args["a3m_file_2"] is not None):
        align_hhblits_output([args["sequence_file_1"],args["sequence_file_2"]], [args["a3m_file_1"],args["a3m_file_2"]], aln_res_file, info_res_file)
    elif (args["sequence_file_1"] is not None) and (args["sequence_file_2"] is not None):
        align_one_hot([args["sequence_file_1"],args["sequence_file_2"]], aln_res_file, info_res_file)
