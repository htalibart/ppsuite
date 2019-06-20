import argparse

from call_compotts import *
from manage_positions import *
from align_msas import *
import files_management as fm
import time


# TODO structure de fichiers


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', '--potts_model_1', help="Potts model 1")
    parser.add_argument('-p2', '--potts_model_2', help="Potts model 2")
    parser.add_argument('-s1', '--sequence_file_1', help="Sequence file 1")
    parser.add_argument('-s2', '--sequence_file_2', help="Sequence file 2")
    parser.add_argument('-h1', '--a3m_file_1', help="HH-blits output file 1")
    parser.add_argument('-h2', '--a3m_file_2', help="HH-blits output file 2")
    parser.add_argument('-o', '--output_folder', help="Output folder")
    parser.add_argument('-of', '--align_train_msas', help="Align MSAs that were used to train the Potts Model using positions aligned by ComPotts", action='store_true')
    parser.add_argument('-os', '--disp_aligned_sequences', help="Display aligned sequences using positions aligned by ComPotts", action='store_true')
    parser.add_argument('-r', '--rescaling_function', help="Rescaling function for Potts model parameters.", default="identity", choices=('identity', 'original_rescaling'))
    parser.add_argument('-nw', '--no_w', help="Don't use w scores", action='store_true')
    parser.add_argument('-m', '--mode', help="Mode", choices=('msgpack', 'hhblits', 'one_hot', 'one_seq_ccmpred'), default='one_seq_ccmpred')
    args = vars(parser.parse_args())


    output_folder = args["output_folder"]
    if output_folder is None:
        general_output_folder = fm.create_folder("output_compotts")
        output_folder = general_output_folder+time.strftime("%Y%m%d-%H%M%S")
    if output_folder[-1]!='/':
        output_folder+='/'
    fm.create_folder(output_folder)

    no_kwargs = ["potts_model_1", "potts_model_2", "sequence_file_1", "sequence_file_2", "a3m_file_1", "a3m_file_2", "output_folder", "mode", "no_w"]
    arguments = {}
    for key in args.keys():
        if key not in no_kwargs:
            arguments[key]=args[key]

    if args['no_w']:
        arguments["w_threshold"]=float('inf')
        arguments["use_w"]=False
    else:
        arguments["use_w"]=True


    if args['mode']=='msgpack': # alignement de deux fichiers msgpack seulement
        if (args["potts_model_1"] is not None) and (args["potts_model_2"] is not None):
            align_two_potts_models_from_files([args["potts_model_1"], args["potts_model_2"]], output_folder, **arguments)
        else:
            print("Need msgpack files")

    else: # tout le reste (fichiers fasta, ...)
        seq_files = [args["sequence_file_1"], args["sequence_file_2"]]

        # récupération des objects ComPotts
        if args['mode']=='hhblits':
            if (args["sequence_file_1"] is not None) and (args["sequence_file_2"] is not None) and (args["a3m_file_1"] is not None) and (args["a3m_file_2"] is not None):
                compotts_objects = [ComPotts_Object.from_hhblits_output(sf, a3m, output_folder, **arguments) for sf, a3m in zip(seq_files, [args["a3m_file_1"], args["a3m_file_2"]])]
            else:
                print("Need sequence files and a3m files")

        elif args['mode']=='one_hot':
            if (args["sequence_file_1"] is not None) and (args["sequence_file_2"] is not None):
                compotts_objects = [ComPotts_Object.from_seq_file_to_one_hot(sf, output_folder, **arguments) for sf in seq_files]
            else:
                print("Need sequence files")


        elif args['mode']=='one_seq_ccmpred':
            if (args["sequence_file_1"] is not None) and (args["sequence_file_2"] is not None):
                compotts_objects = [ComPotts_Object.from_seq_file_via_ccmpred(sf, output_folder, **arguments) for sf in seq_files]
            else:
                print("Need sequence files")

 
        # alignement
        aligned_positions, infos_solver = align_two_objects(compotts_objects, output_folder, **arguments)


        # on fait des trucs avec les positions alignées
        if args["align_train_msas"]:
            output_msa = output_folder+'_'.join(o.name for o in compotts_objects)+".fasta"
            get_msas_aligned(aligned_positions, [o.train_msa for o in compotts_objects], output_msa)

        if args["disp_aligned_sequences"]:
            seqs_aligned = get_seqs_aligned(aligned_positions, compotts_objects)
            for s in seqs_aligned:
                print(s)
