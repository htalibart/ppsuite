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
        sub_output_folder = create_folder(output_folder+obj1.name+"_"+obj2.name+"/")
        if not os.path.isfile(fm.get_aln_res_filename(sub_output_folder)):
            aligned_positions, info_solver = align_two_objects([obj1, obj2], sub_output_folder, **params)
        return Compotts_Object.from_merge(obj1, obj2, aligned_positions, **kwargs)


def multiple_alignment(seq_files, a3m_files, output_folder, **kwargs):
    objs_dict = collections.OrderedDict()
    for seq_file, a3m_file in zip(seq_files, a3m_files):
        obj = ComPotts_Object.from_hhblits_output(seq_file, a3m_file, **kwargs)
        objs_dict[obj.name] = obj
    first_node = cluster_compotts_objects(objs_dict)
    get_subalignment(first_node, obj_dict, output_folder)
