import argparse
import sys
import time
import pathlib
import subprocess

from compotts.call_compotts import *
from compotts.manage_positions import *
from makemsa.__main__ import *
from vizpotts.vizpotts import *
import comutils.files_management as fm


def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-p1', '--potts_model_file_1', help="Potts model 1", type=pathlib.Path)
    parser.add_argument('-p2', '--potts_model_file_2', help="Potts model 2", type=pathlib.Path)
    parser.add_argument('-f1', '--feature_folder_1', help="Folder containing files for sequence 1", type=pathlib.Path)
    parser.add_argument('-f2', '--feature_folder_2', help="Folder containing files for sequence 2", type=pathlib.Path)
    parser.add_argument('-gf1', '--guess_folder_1', help=argparse.SUPPRESS, type=pathlib.Path)
    parser.add_argument('-gf2', '--guess_folder_2', help=argparse.SUPPRESS, type=pathlib.Path)
    parser.add_argument('-o', '--output_folder', help="Output folder (if not specified : output_compotts/[DATE]/)", type=pathlib.Path)

    # options alignement
    parser.add_argument('-nw', '--no_w', help="Don't use w scores (default : False)", action='store_true')
    parser.add_argument('-nv', '--no_v', help="Don't use v scores (default : False)", action='store_true')
    parser.add_argument('-wt', '--w_threshold_method', help="w threshold method. Couplings that have a Frobenius norm below the threshold are not considered by ComPotts", default="no_threshold") # TODO checker si c'est bien fait avant le rescaling
    parser.add_argument('--v_rescaling_function', help="Rescaling function for the v parameters of the Potts model. (default : no rescaling (identity))", default="identity")
    parser.add_argument('--w_rescaling_function', help="Rescaling function for the w parameters of the Potts model. (default : no rescaling (identity))", default="identity")
    parser.add_argument('--v_shift', help=argparse.SUPPRESS, type=float, default=3)
    parser.add_argument('--wijab_threshold', help="If |wijab|<wijab_threshold, wijab is set to 0 if using rescaling function threshold_on_wijab", type=float, default=0)
    parser.add_argument('--v_rescaling_tau', help="Tau parameter for rescaling function simulate_uniform_pc_on_v", type=float, default=0.5)
    parser.add_argument('--no_v_back_to_scale', help="Don't put v back to old norm after simulate_uniform_pc_on_v", default=False, action='store_true')
    parser.add_argument('--w_rescaling_tau', help="Tau parameter for rescaling function simulate_uniform_pc_on_w", type=float, default=0.5)
    parser.add_argument('--add_pseudo_w', help=argparse.SUPPRESS, action='store_true')
    parser.add_argument('--rescale_wij', help=argparse.SUPPRESS, action='store_true')
    parser.add_argument('--w_submat_tau', help=argparse.SUPPRESS, type=float, default=0.05)
    parser.add_argument('--beta_softmax_w', help="Beta rescaling parameter to simulate uniform pseudo-counts on w through softmax", type=float, default=10)
    # parser.add_argument('-vwc', '--vw_coeff_method', help=argparse.SUPPRESS, default="arbitrary_1_1") # v w coeff method
    parser.add_argument('--alpha_w', help="coefficient before w score", default=1, type=float)
    #parser.add_argument('-gc', '--gap_cost_method', help=argparse.SUPPRESS, default="arbitrary_8_0") # gap costs method
    parser.add_argument('-go', '--gap_open', help="gap open", default=8, type=float) # gap costs method
    parser.add_argument('-ge', '--gap_extend', help="gap extend", default=0, type=float) # gap costs method
    #parser.add_argument('--gap_auto_coeff', help=argparse.SUPPRESS, type=float, default=2.5) # gap costs method (soon to be deprecated)
    # solver options
    parser.add_argument('-t', '--t_limit', help="solver : time limit in seconds (default : 36000)", type=float, default=36000)
    parser.add_argument('-lit', '--iter_limit_param', help="solver : nb Lagrangian iterations (default : 1000000000)", type=int, default=1000000000)
    #parser.add_argument('-e', '--epsilon', help="solver : precision", type=float, default=1)
    parser.add_argument('-e', '--precision_method', help="solver : precision method (default : similarity_0.005)", default="similarity_0.005") # TODO documenter toutes les précisions
    parser.add_argument('-ga', '--gamma', help="solver : gamma (default : 1.0)", type=float, default=1.0)
    parser.add_argument('-th', '--theta', help="solver : theta (default : 0.9)", type=float, default=0.9)
    parser.add_argument('-stpz', '--stepsize_min', help="solver : stepsize_min (default : 0.000000005)", type=float, default=0.000000005)
    parser.add_argument('-stpm', '--nb_non_increasing_steps_max', help="solver : nb_non_increasing_steps_max (default : 500)", type=int, default=500)
    parser.add_argument('-sim_min', '--sim_min', help='if similarity score is below sim_min, solver considers that the Potts models are not similar and stops. (default : 0)', type=float, default=0)

    # autres options
    parser.add_argument('-ali', '--call_aliview', help=argparse.SUPPRESS, action='store_true')
    parser.add_argument('-oaln', '--get_training_sets_fasta_aln', help="Get training sets alignment in a fasta file", action='store_true')
    parser.add_argument('-osaln', '--get_sequences_fasta_aln', help="Get sequences alignment in a fasta file", action='store_true')

    args = vars(parser.parse_args(args))

    if args["gap_extend"]>args["gap_open"]:
        raise Exception("Gap extend must be smaller than gap open (convergence issues with the solver)")


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


    args['v_back_to_scale'] = not args['no_v_back_to_scale']

    # INSTANCIATING OBJECTS
    compotts_objects = []
    for k in range(1,3):
        if args["potts_model_file_"+str(k)] is not None:
            feature_folder = pathlib.Path(tempfile.mkdtemp())
            shutil.copy(args["potts_model_file_"+str(k)], feature_folder/"potts_model.mrf")
            obj = Potts_Object.from_folder(feature_folder, **args)
            compotts_objects.append(obj)
        elif args["feature_folder_"+str(k)] is not None:
            obj = Potts_Object.from_folder(args["feature_folder_"+str(k)], **args)
            compotts_objects.append(obj)
        elif args["guess_folder_"+str(k)] is not None:
            obj = Potts_Object.guess_from_folder(args["guess_folder_"+str(k)], **args)
            compotts_objects.append(obj)
        else:
            raise Exception("Need input "+str(k))

    
    # WRITE README
    fm.write_readme(output_folder, **args)


    # ALIGNMENT
    aligned_positions, infos_solver = align_two_objects(compotts_objects, output_folder, **args)
    print("Total time : "+str(infos_solver["total_compotts_time"]))


    # MAYBE DO SOMETHING WITH THE ALIGNMENT

    if len(aligned_positions)>0:

        # GIVE ALIGNED POSITIONS FOR THE ORIGINAL ALIGNMENTS
        if all((o.mrf_pos_to_aln_pos is not None) for o in compotts_objects):
            original_positions = get_initial_positions(aligned_positions, {"pos_ref":compotts_objects[0].mrf_pos_to_aln_pos, "pos_2":compotts_objects[1].mrf_pos_to_aln_pos})
            fm.write_positions_to_csv(original_positions, output_folder/("aln_original.csv"))

        # GIVE ALIGNED POSITIONS FOR THE SEQUENCES
        if all((o.mrf_pos_to_seq_pos is not None) for o in compotts_objects):
            sequence_positions = get_initial_positions(aligned_positions, {"pos_ref":compotts_objects[0].mrf_pos_to_seq_pos, "pos_2":compotts_objects[1].mrf_pos_to_seq_pos})
            fm.write_positions_to_csv(sequence_positions, output_folder/("aln_sequences.csv"))
 
        if all((o.sequence is not None) for o in compotts_objects) and args["get_sequences_fasta_aln"]:
            output_fasta_file = output_folder/("aligned_sequences.fasta")
            get_seqs_aligned_in_fasta_file(aligned_positions, compotts_objects, output_fasta_file)


        # ALIGN TRAINING MSAS
        if all((o.aln_train is not None) for o in compotts_objects) and args["get_training_sets_fasta_aln"]:
            output_msa = output_folder/("aligned_training_sets.fasta")
            get_msas_aligned(aligned_positions, [o.aln_train for o in compotts_objects], output_msa)
            if args["call_aliview"]:
                cmd = "aliview "+str(output_msa)
                subprocess.Popen(cmd, shell=True).wait()

           
        return {"compotts_objects": compotts_objects, "aligned_positions":aligned_positions, "infos_solver":infos_solver}



if __name__=="__main__":
    main()
