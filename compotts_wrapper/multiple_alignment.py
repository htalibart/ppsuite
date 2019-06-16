from scipy.cluster.hierarchy import *
from matplotlib import pyplot as plt
import collections
from compotts_wrapper.compotts_object import *

# TODO tester

def distance_metric(obj1, obj2):
    return 1-seq_identity(obj1.real_seq, obj2.real_seq)


def cluster_compotts_objects(objs, draw_dendro=True):
    Z = linkage([distance_metric(compotts_objects[seq_names[j]], compotts_objects[seq_names[k]]) for j in range(len(compotts_objects)) for k in range(j+1, len(compotts_objects))])
    if show_figure:
        plt.figure()
        dendro = dendrogram(Z, labels=seq_names, leaf_font_size=8, distance_sort='ascending')
        plt.show()
    return to_tree(Z)


def get_subalignment(node, compotts_objects, output_folder):
    if node.is_leaf():
        return compotts_objects.keys()[node.id]
    else:
        obj1 = get_subalignment(node.left, compotts_objects, output_folder)
        obj2 = get_subalignment(node.right, compotts_objects, output_folder)
        res_aln_file = output_folder+obj1.name+"_"+obj2.name+"_aln.csv"
        res_info_file = output_folder+obj1.name+"_"+obj2.name+"_info.csv"
        if not os.path.isfile(res_aln_file):
            align_two_objects([obj1, obj2], res_aln_file, res_info_file, **params)
        return Compotts_Object.from_merge(obj1, obj2, res_aln_file, **kwargs)


def multiple_alignment(seq_files, a3m_files, output_folder, **kwargs):
    objs_dict = collections.OrderedDict()
    for seq_file in seq_files:
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, **kwargs)
        objs_dict[obj.name] = obj
    first_node = cluster_compotts_objects(objs_dict)
    get_subalignment(first_node, obj_dict, output_folder)
