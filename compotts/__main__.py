import argparse

from compotts.call_compotts import *
from compotts.manage_positions import *
from compotts.align_msas import *
import basic_modules.files_management as fm
import time

import json

# TODO structure de fichiers

def write_readme(folder, **kwargs):
    with open(folder+"README.txt", 'w') as f:
        json.dump(kwargs, f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', '--potts_model_1', help="Potts model 1")
    parser.add_argument('-p2', '--potts_model_2', help="Potts model 2")
    parser.add_argument('-s1', '--sequence_file_1', help="Sequence file 1")
    parser.add_argument('-s2', '--sequence_file_2', help="Sequence file 2")
    parser.add_argument('-h1', '--a3m_file_1', help="HH-blits output file 1")
    parser.add_argument('-h2', '--a3m_file_2', help="HH-blits output file 2")
    parser.add_argument('-f1', '--input_folder_1', help="Folder containing files for sequence 1")
    parser.add_argument('-f2', '--input_folder_2', help="Folder containing files for sequence 2")
    parser.add_argument('-o', '--output_folder', help="Output folder")
    parser.add_argument('-of', '--align_train_msas', help="Align MSAs that were used to train the Potts Model using positions aligned by ComPotts", action='store_true') # devrait être fait par défaut...?
    parser.add_argument('-os', '--align_sequences', help="Display aligned sequences using positions aligned by ComPotts", action='store_true')
    parser.add_argument('-r', '--rescaling_function', help="Rescaling function for Potts model parameters.", default="identity", choices=('identity', 'original_rescaling'))
    parser.add_argument('-nw', '--no_w', help="Don't use w scores", action='store_true')
    parser.add_argument('-wt', '--w_threshold_method', help="w threshold method. Couplings that have a Frobenius norm below the threshold are not considered by ComPotts", default="percentile_80") # TODO checker si c'est bien fait avant le rescaling
    parser.add_argument('-m', '--mode', help="Mode", choices=('msgpack', 'hhblits', 'one_hot', 'one_seq_ccmpred'), default='one_seq_ccmpred')
    args = vars(parser.parse_args())


    output_folder = args["output_folder"]
    if output_folder is None:
        general_output_folder = fm.create_folder("output_compotts")
        output_folder = general_output_folder+time.strftime("%Y%m%d-%H%M%S")
    if output_folder[-1]!='/':
        output_folder+='/'
    fm.create_folder(output_folder)

    for k in range(1,3):
        if args["input_folder_"+str(k)] is not None:
            if args["input_folder_"+str(k)][-1]!='/':
                args["input_folder_"+str(k)]+='/'

    no_kwargs = ["potts_model_1", "potts_model_2", "sequence_file_1", "sequence_file_2", "a3m_file_1", "a3m_file_2", "output_folder", "mode", "no_w"] # TODO voir si utile
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
        for k in range(1,3):
            if (args["potts_model_"+str(k)] is None) and (args["input_folder_"+str(k)] is not None):
                args["potts_model_"+str(k)] = fm.get_potts_model_file_from_folder(args["input_folder_"+str(k)])
        # TODO rescaling
        if (args["potts_model_1"] is not None) and (args["potts_model_2"] is not None):
            align_two_potts_models_from_files([args["potts_model_1"], args["potts_model_2"]], output_folder, **arguments)
        else:
            print("Need msgpack files")



    else: # tout le reste (fichiers fasta, ...)
        for k in range(1,3):
            if (args["sequence_file_"+str(k)] is None) and (args["input_folder_"+str(k)] is not None):
                args["sequence_file_"+str(k)] = fm.get_sequence_file_from_folder(args["input_folder_"+str(k)])

        seq_files = [args["sequence_file_1"], args["sequence_file_2"]]

        # récupération des objects ComPotts
        if args['mode']=='hhblits':
            for k in range(1,3):
                if (args["a3m_file_"+str(k)] is None) and (args["input_folder_"+str(k)] is not None):
                    args["a3m_file_"+str(k)] = fm.get_a3m_file_from_folder(args["input_folder_"+str(k)])
                if (args["potts_model_"+str(k)] is None) and (args["input_folder_"+str(k)] is not None):
                    args["potts_model_"+str(k)] = fm.get_potts_model_file_from_folder(args["input_folder_"+str(k)])

            if (args["sequence_file_1"] is not None) and (args["sequence_file_2"] is not None) and (args["a3m_file_1"] is not None) and (args["a3m_file_2"] is not None):
                compotts_objects = []
                for k in range(2):
                    obj = ComPotts_Object.from_hhblits_output(seq_files[k], args["a3m_file_"+str(k+1)], output_folder, mrf_file=args["potts_model_"+str(k+1)], **arguments)
                    compotts_objects.append(obj)
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


        write_readme(output_folder, **arguments)

 
        # alignement
        aligned_positions, infos_solver = align_two_objects(compotts_objects, output_folder, **arguments)

        print("Total time : "+str(infos_solver["total_compotts_time"]))

        # on fait des trucs avec les positions alignées
        if args["align_train_msas"]:
            output_msa = output_folder+'_'.join(o.name for o in compotts_objects)+"_train_msas.fasta"
            get_msas_aligned(aligned_positions, [o.train_msa for o in compotts_objects], output_msa)

        if args["align_sequences"]:
            output_fasta_file = output_folder+'_'.join(o.name for o in compotts_objects)+"_aligned_sequences.fasta"
            get_seqs_aligned_in_fasta_file(aligned_positions, compotts_objects, output_fasta_file)


if __name__=="__main__":
    main()
