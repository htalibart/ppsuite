import argparse
import sys
import time
import pathlib
import subprocess

from ppalign.call_ppalign import *
from ppalign.manage_positions import *
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
    parser.add_argument('-o', '--output_folder', help="Output folder (if not specified : output_ppalign/[DATE]/)", type=pathlib.Path)


    # options to pre-process potts models
    parser.add_argument('--w_percent', help="% couplings considered (wij with lowest norms are set to 0)", default=100, type=float)
    parser.add_argument('--v_rescaling_function', help=argparse.SUPPRESS, default="simulate_uniform_pc_on_v")
    parser.add_argument('--v_rescaling_tau', help="Tau parameter for field smoothing", type=float, default=0.4)
    parser.add_argument('--v_back_to_scale', help=argparse.SUPPRESS, default=False, action='store_true')
    parser.add_argument('--w_rescaling_function', help=argparse.SUPPRESS, default="simulate_uniform_pc_on_w")
    parser.add_argument('--w_rescaling_tau', help="Tau parameter for coupling smoothing", type=float, default=0.4)
    parser.add_argument('--beta_softmax_w', help="Softmax base parameter for coupling smoothing", type=float, default=10)
    parser.add_argument('--w_back_to_scale', help=argparse.SUPPRESS, default=False, action='store_true')
    parser.add_argument('--insert_null_at_trimmed', help="insert null parameters at trimmed positions", default=False, action='store_true')

    # options alignement
    parser.add_argument('-nw', '--no_w', help="Don't use w scores (default : False)", action='store_true')
    parser.add_argument('-nv', '--no_v', help="Don't use v scores (default : False)", action='store_true')
    parser.add_argument('--alpha_w', help="coefficient before w score", default=1, type=float)
    parser.add_argument('--offset_v', help="score offset for v parameters", default=0, type=float)
    parser.add_argument('--exp', help="scalar product of the exp instead of simple scalar product", action='store_true', default=False)
    parser.add_argument('--remove_v0', help="remove background v0", action='store_true', default=False)
    parser.add_argument('--rescale_removed_v0', help="rescale background v0", action='store_true', default=False)
    parser.add_argument('-go', '--gap_open', help="gap open", default=14, type=float)
    parser.add_argument('-ge', '--gap_extend', help="gap extend", default=0, type=float)

    # solver options
    parser.add_argument('-t', '--t_limit', help="solver : time limit in seconds (default : 36000)", type=float, default=36000)
    parser.add_argument('-lit', '--iter_limit_param', help="solver : nb Lagrangian iterations (default : 1000000000)", type=int, default=1000000000)
    parser.add_argument('--epsilon_sim', help="solver : max 2*(UB-LB)/(s(A,A)+s(B,B)) (default : 0.005)", default=0.005, type=float)
    parser.add_argument('-ga', '--gamma', help="solver : gamma (default : 1.0)", type=float, default=1.0)
    parser.add_argument('-th', '--theta', help="solver : theta (default : 0.9)", type=float, default=0.9)
    parser.add_argument('-stpz', '--stepsize_min', help="solver : stepsize_min (default : 0.000000005)", type=float, default=0.000000005)
    parser.add_argument('-stpm', '--nb_non_increasing_steps_max', help="solver : nb_non_increasing_steps_max (default : 500)", type=int, default=500)
    parser.add_argument('-sim_min', '--sim_min', help='if similarity score is below sim_min, solver considers that the Potts models are not similar and stops. (default : -100)', type=float, default=-100)

    # other output options
    parser.add_argument('-ali', '--call_aliview', help=argparse.SUPPRESS, action='store_true')
    parser.add_argument('-oaln', '--get_training_sets_fasta_aln', help="Get training sets alignment in a fasta file", action='store_true')
    parser.add_argument('-osaln', '--get_sequences_fasta_aln', help="Get sequences alignment in a fasta file", action='store_true')

    args = vars(parser.parse_args(args))

    if args["gap_extend"]>args["gap_open"]:
        raise Exception("Gap extend must be smaller than gap open (convergence issues with the solver)")


    # CREATE_FOLDER IF NOT EXISTING
    if args["output_folder"] is None:
        general_output_folder = fm.create_folder("output_ppalign")
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


    # INSTANCIATING OBJECTS
    objects = []
    temp_folders = []
    for k in range(1,3):
        if args["potts_model_file_"+str(k)] is not None:
            feature_folder = pathlib.Path(tempfile.mkdtemp())
            shutil.copy(str(args["potts_model_file_"+str(k)]), str(feature_folder/"potts_model.mrf"))
            obj = Potts_Object.from_folder(feature_folder, **args)
            objects.append(obj)
            temp_folders.append(feature_folder)
        elif args["feature_folder_"+str(k)] is not None:
            obj = Potts_Object.from_folder(args["feature_folder_"+str(k)], **args)
            objects.append(obj)
        elif args["guess_folder_"+str(k)] is not None:
            obj = Potts_Object.guess_from_folder(args["guess_folder_"+str(k)], **args)
            objects.append(obj)
        else:
            raise Exception("Need input "+str(k))

    # EXP IF NEEDED
    for obj in objects:
        if args["exp"]:
            obj.potts_model = get_rescaled_potts_model(obj.potts_model, "exponential", "exponential", args["use_w"])
    
    
    # INSERT NULL AT TRIMMED
    for obj in objects:
        if args["insert_null_at_trimmed"]:
            if (args["feature_folder_1"] is None) or (args["feature_folder_2"]) is None:
                raise Exception("folders must be specified to insert null columns")
            obj.insert_null_at_trimmed()

    # WRITE README
    fm.write_readme(output_folder, **args)


    # ALIGNMENT
    aligned_positions, infos_solver = align_two_objects(objects, output_folder, **args)
    print("Total time : "+str(infos_solver["total_time"]))


    # MAYBE DO SOMETHING WITH THE ALIGNMENT

    if len(aligned_positions)>0:

      #  # GIVE ALIGNED POSITIONS FOR THE ORIGINAL ALIGNMENTS
      #  if all((o.mrf_pos_to_aln_pos is not None) for o in objects):
      #      original_positions = get_initial_positions(aligned_positions, {"pos_ref":objects[0].mrf_pos_to_aln_pos, "pos_2":objects[1].mrf_pos_to_aln_pos})
      #      fm.write_positions_to_csv(original_positions, output_folder/("aln_original.csv"))

        # GIVE ALIGNED POSITIONS FOR THE SEQUENCES
        if all((o.mrf_pos_to_seq_pos is not None) for o in objects):
            sequence_positions = get_initial_positions(aligned_positions, {"pos_ref":objects[0].mrf_pos_to_seq_pos, "pos_2":objects[1].mrf_pos_to_seq_pos})
            fm.write_positions_to_csv(sequence_positions, output_folder/("aln_sequences.csv"))
 
        if all((o.sequence is not None) for o in objects) and args["get_sequences_fasta_aln"]:
            output_fasta_file = output_folder/("aligned_sequences.fasta")
            get_seqs_aligned_in_fasta_file(aligned_positions, objects, output_fasta_file)


        # ALIGN TRAINING MSAS
        if all((o.aln_train is not None) for o in objects) and args["get_training_sets_fasta_aln"]:
            output_msa = output_folder/("aligned_training_sets.fasta")
            get_msas_aligned(aligned_positions, [o.aln_train for o in objects], output_msa)
            if args["call_aliview"]:
                cmd = "aliview "+str(output_msa)
                subprocess.Popen(cmd, shell=True).wait()

        # REMOVE TEMPORARY FOLDERS
        for temp_folder in temp_folders: 
            if temp_folder.is_dir():
                shutil.rmtree(str(temp_folder))

        return {"objects": objects, "aligned_positions":aligned_positions, "infos_solver":infos_solver}



if __name__=="__main__":
    main()
