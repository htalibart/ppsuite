import argparse
import pathlib
import os
import json
from scipy.cluster.hierarchy import *
from matplotlib import pyplot as plt
import collections
from compotts.compotts_object import *
from compotts.call_compotts import *
import basic_modules.files_management as fm

# TODO README

def seq_id_distance(obj1, obj2):
    return 1-seq_identity(obj1.real_seq, obj2.real_seq)

def inverse_seq_id_distance(obj1, obj2):
    return seq_identity(obj1.real_seq, obj2.real_seq)


def cluster_compotts_objects(objs, distance_metric, draw_dendro=False):
    seq_names = list(objs.keys())
    Z = linkage([eval(distance_metric)(objs[seq_names[j]], objs[seq_names[k]]) for j in range(len(objs)) for k in range(j+1, len(objs))])
    if draw_dendro:
        plt.figure()
        dendro = dendrogram(Z, labels=seq_names, leaf_font_size=8, distance_sort='ascending')
        plt.show()
    return to_tree(Z)


def get_subalignment(node, compotts_objects, output_folder, **kwargs):
    if node.is_leaf():
        return compotts_objects[list(compotts_objects.keys())[node.id]]
    else:
        obj1 = get_subalignment(node.left, compotts_objects, output_folder, **kwargs)
        obj2 = get_subalignment(node.right, compotts_objects, output_folder, **kwargs)
        sub_output_folder = fm.create_folder(os.path.join(output_folder,obj1.name+"_"+obj2.name))
        if not os.path.isfile(fm.get_aln_res_file_name(sub_output_folder)):
            aligned_positions, info_solver = align_two_objects([obj1, obj2], sub_output_folder, **kwargs)
        else:
            aligned_positions = fm.get_aligned_positions_dict_from_compotts_output_file(fm.get_aln_res_file_name(sub_output_folder))
        return ComPotts_Object.from_merge(obj1, obj2, aligned_positions, output_folder, **kwargs)


def multiple_alignment(protein_folders, output_folder, distance_metric="seq_id_distance", **kwargs):
    objs_dict = collections.OrderedDict()
    for pf in protein_folders:
        seq_file = fm.get_sequence_file_from_folder(pf)
        a3m_file = fm.get_a3m_file_from_folder(pf)
        mrf_file = fm.get_potts_model_file_from_folder(pf)
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, output_folder, input_folder=pf, mrf_file=mrf_file, **kwargs)
        objs_dict[obj.name] = obj
    first_node = cluster_compotts_objects(objs_dict, distance_metric=distance_metric)
    get_subalignment(first_node, objs_dict, output_folder, **kwargs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folders', help="folders for each protein (with sequence files and a3m files)", nargs='+', default=[])
    parser.add_argument('-of', '--output_folder', help="output folder", default="output_multiple", type=pathlib.Path)
    parser.add_argument('-r', '--rescaling_function', help="Rescaling function for Potts model parameters.", default="identity", choices=('identity', 'original_rescaling', 'symmetric_relu_like'))
    parser.add_argument('-nw', '--no_w', help="Don't use w scores", action='store_true')
    parser.add_argument('-wt', '--w_threshold_method', help="w threshold method. Couplings that have a Frobenius norm below the threshold are not considered by ComPotts", default="no_threshold") # TODO checker si c'est bien fait avant le rescaling
    parser.add_argument('-go', '--gap_open', help="gap open", type=float, default=0)
    parser.add_argument('-ge', '--gap_extend', help="gap extend", type=float, default=0)
    parser.add_argument('-dm', '--distance_metric', help="distance metric", choices=("seq_id_distance", "inverse_seq_id_distance"), default="seq_id_distance")
    args = vars(parser.parse_args())

    if not os.path.isdir(args["output_folder"]):
        os.mkdir(args["output_folder"])

    fm.write_readme(args["output_folder"], **args)

    multiple_alignment(args["folders"], args["output_folder"], distance_metric=args["distance_metric"],
    use_w=(not args["no_w"]), w_threshold_method=args["w_threshold_method"], gap_open=args["gap_open"], gap_extend=args["gap_extend"])


if __name__=="__main__":
    main()
