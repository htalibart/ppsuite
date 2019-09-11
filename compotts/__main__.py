import argparse
import sys
import time
import pathlib

from compotts.call_compotts import *
from compotts.manage_positions import *
from compotts.align_msas import *
from vizpotts.vizpotts import *
import basic_modules.files_management as fm


def handle_args_for_obj(args, k):
    """ Returns a dictionary containing arguments @args except only those that are relevant for object @k and cleaning the names """
    new_args = {}
    for key in args:
        if str(key).endswith(str(k+1)):
            new_key = key[:-len("_"+str(k))]
            new_args[new_key] = args[key]
        elif not key.endswith(str((k+1)%2+1)):
            new_args[key] = args[key]
    return new_args




def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', '--potts_model_file_1', help="Potts model 1", type=pathlib.Path)
    parser.add_argument('-p2', '--potts_model_file_2', help="Potts model 2", type=pathlib.Path)
    parser.add_argument('-s1', '--sequence_file_1', help="Sequence file 1", type=pathlib.Path)
    parser.add_argument('-s2', '--sequence_file_2', help="Sequence file 2", type=pathlib.Path)
    parser.add_argument('-h1', '--a3m_file_1', help="HH-blits output file 1", type=pathlib.Path)
    parser.add_argument('-h2', '--a3m_file_2', help="HH-blits output file 2", type=pathlib.Path)
    parser.add_argument('-f1', '--input_folder_1', help="Folder containing files for sequence 1", type=pathlib.Path)
    parser.add_argument('-f2', '--input_folder_2', help="Folder containing files for sequence 2", type=pathlib.Path)
    parser.add_argument('-aln1', '--aln_fasta_1', help="Alignment file in fasta format 1", type=pathlib.Path)
    parser.add_argument('-aln2', '--aln_fasta_2', help="Alignment file in fasta format 2", type=pathlib.Path)
    parser.add_argument('-o', '--output_folder', help="Output folder", type=pathlib.Path)
    parser.add_argument('-r', '--rescaling_function', help="Rescaling function for Potts model parameters.", default="identity", choices=('identity', 'original_rescaling', 'symmetric_relu_like', 'shifted_relu'))
    parser.add_argument('-n', '--nb_sequences', help="Number of sequences in the MRF training alignment", default=1000, type=int)
    parser.add_argument('-m1', '--mrf_type_1', help="Mode 1", choices=('standard', 'one_hot', 'one_submat'), default='standard')
    parser.add_argument('-m2', '--mrf_type_2', help="Mode 2", choices=('standard', 'one_hot', 'one_submat'), default='standard')
    parser.add_argument('-nm1', '--name_1', help="Name for alignee 1")
    parser.add_argument('-nm2', '--name_2', help="Name for alignee 2")

    # trimal
    parser.add_argument('-trimgt', '--trimal_gt', help="trimal gt", default=0.8, type=float)
    parser.add_argument('-trimcons', '--trimal_cons', help="trimal cons", default=60, type=float)

    # options alignement
    parser.add_argument('-nw', '--no_w', help="Don't use w scores", action='store_true')
    parser.add_argument('-nv', '--no_v', help="Don't use v scores", action='store_true')
    parser.add_argument('-wt', '--w_threshold_method', help="w threshold method. Couplings that have a Frobenius norm below the threshold are not considered by ComPotts", default="no_threshold") # TODO checker si c'est bien fait avant le rescaling
    parser.add_argument('-vwc', '--vw_coeff_method', help="vw coeff method", default="arbitrary_1_1") # TODO doc
    parser.add_argument('-gc', '--gap_cost_method', help="gap costs method", default="arbitrary_8_0") # TODO doc

    # solver options
    parser.add_argument('-t', '--t_limit', help="solver : time limit", type=float, default=36000)
    parser.add_argument('-lit', '--iter_limit_param', help="solver : nb Lagrangian iterations", type=int, default=1000000000)
    #parser.add_argument('-e', '--epsilon', help="solver : precision", type=float, default=1)
    parser.add_argument('-e', '--epsilon_method', help="solver : precision method", default="arbitrary_1")
    parser.add_argument('-ga', '--gamma', help="solver : gamma", type=float, default=1.0)
    parser.add_argument('-th', '--theta', help="solver : theta", type=float, default=0.9)
    parser.add_argument('-stpz', '--stepsize_min', help="solver : stepsize_min", type=float, default=0.000000005)
    parser.add_argument('-stpm', '--nb_non_increasing_steps_max', help="solver : nb_non_increasing_steps_max", type=int, default=500)
    

    # autres options
    parser.add_argument('-ali', '--call_aliview', help="Call aliview at the end", action='store_true')


    # CCMpredPy options
    parser.add_argument('--pc_count', help="CCMpred : Specify number of pseudocounts (default : 1000)", default=1000)
    parser.add_argument('--reg_lambda_pair_factor', help="CCMpred : Regularization parameter for pair potentials (L2 regularization with lambda_pair = lambda_pair-factor * scaling) [CCMpred default: 0.2, our default : 30]", default=30)

    # options related to CCMpredPy
    #parser.add_argument('--pc_count_factor', help="Number of pseudocounts for CCMpredPy will be (pc_count_factor)*(nb_sequences). Default : pc_count_factor=1", type=int, default=1)


    # options related to HHblits
    parser.add_argument('-d', '--hhblits_database', help="Database for HHblits", type=pathlib.Path)


    args = vars(parser.parse_args(args))


    # CREATE_FOLDER IF NOT EXISTING
    if args["output_folder"] is None:
        general_output_folder = fm.create_folder("output_compotts")
        output_folder = general_output_folder / time.strftime("%Y%m%d-%H%M%S")
    else:
        output_folder = args["output_folder"]
    fm.create_folder(output_folder)
    del args["output_folder"]


    # HANDLING NO W / NO V ARGUMENTS
    if args["no_w"]:
        args["w_threshold"]=float('inf')
        args["use_w"]=False
    else:
        args["use_w"]=True
    del args["no_w"]

    args["use_v"]= not args["no_v"]
    del args["no_v"]


    # CREATING COMPOTTS OBJECTS
    compotts_objects = [] 
    for k in range(2):
        new_args = handle_args_for_obj(args, k)
        obj = ComPotts_Object(**new_args)
        compotts_objects.append(obj) 


    # WRITE README
    fm.write_readme(output_folder, **args)


    # ALIGNMENT
    aligned_positions, infos_solver = align_two_objects(compotts_objects, output_folder, **args)
    print("Total time : "+str(infos_solver["total_compotts_time"]))



    # MAYBE DO SOMETHING WITH THE ALIGNMENT

    # ALIGN TRAINING MSAS
    if all((o.training_set is not None) for o in compotts_objects):
        output_msa = output_folder/('_'.join(o.name for o in compotts_objects)+"_aligned_training_sets.fasta")
        get_msas_aligned(aligned_positions, [o.training_set for o in compotts_objects], output_msa)
        if args["call_aliview"]:
            os.system("aliview "+str(output_msa))
        
    # ALIGN SEQUENCES
    if all((o.sequence is not None) for o in compotts_objects):
        output_fasta_file = output_folder/('_'.join(o.name for o in compotts_objects)+"_aligned_sequences.fasta")
        get_seqs_aligned_in_fasta_file(aligned_positions, compotts_objects, output_fasta_file)


        return {"compotts_objects": compotts_objects, "aligned_positions":aligned_positions, "infos_solver":infos_solver}



if __name__=="__main__":
    main()
