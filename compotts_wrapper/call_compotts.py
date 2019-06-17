from compotts_wrapper.compute_scores import *
from compotts_wrapper.compotts_object import *
from potts_model import *
import ctypes

COMPOTTS_SOLVER = ctypes.CDLL("./compotts_solver.so")
INFINITY = 10000000 



# TODO stocker les scores en utilisant la symétrie pour prendre moins de mémoire
# TODO aln_res_file et info_res_file comme sortie de ComPotts -> variables

def align_two_potts_models(mrfs, output_folder, n_limit_param=INFINITY, iter_limit_param=INFINITY, t_limit=36000, disp_level=1, epsilon=1, v_score_function=scalar_product, w_score_function=scalar_product, gap_open=0, gap_extend=0, w_threshold=0, **kwargs):

    aln_res_file = output_folder+"aln.csv"
    info_res_file = output_folder+"info.csv"

    v_scores = compute_v_scores(*mrfs, v_score_function)
    c_v_scores = ctypes.c_void_p(v_scores.ctypes.data)

    edges_maps = [get_edges_map(mrf, w_threshold) for mrf in mrfs]
    c_int_p = ctypes.POINTER(ctypes.c_int)
    c_edges_mapA = np.ascontiguousarray(edges_maps[0].flatten(), dtype=np.int32).ctypes.data_as(c_int_p)
    c_edges_mapB = np.ascontiguousarray(edges_maps[1].flatten(), dtype=np.int32).ctypes.data_as(c_int_p)

    w_scores = compute_w_scores(*mrfs, *edges_maps, w_score_function)
    c_w_scores = ctypes.c_void_p(w_scores.ctypes.data)

    selfcompA = compute_selfscore(mrfs[0], edges_maps[0], v_score_function, w_score_function)
    selfcompB = compute_selfscore(mrfs[1], edges_maps[1], v_score_function, w_score_function)
    
    COMPOTTS_SOLVER.call_from_python(c_v_scores, c_w_scores, ctypes.c_int(mrfs[0].ncol), ctypes.c_int(mrfs[1].ncol), c_edges_mapA, c_edges_mapB, ctypes.c_double(selfcompA), ctypes.c_double(selfcompB), ctypes.c_double(gap_open), ctypes.c_double(gap_extend), ctypes.c_char_p(aln_res_file.encode('utf-8')), ctypes.c_char_p(info_res_file.encode('utf-8')), ctypes.c_int(iter_limit_param), ctypes.c_int(t_limit), ctypes.c_int(disp_level), ctypes.c_double(epsilon))


def align_two_potts_models_from_files(mrf_files, output_folder, **kwargs):
    mrfs = [Potts_Model.from_msgpack(mrf_file) for mrf_file in mrf_files]
    align_two_potts_models(mrfs, output_folder, **kwargs)


def align_two_objects(objects, output_folder, **kwargs):
    align_two_potts_models([o.mrf for o in objects], output_folder, **kwargs)


def align_hhblits_output(seq_files, a3m_files, output_folder, **kwargs):
    objects = []
    for seq_file, a3m_file in zip(seq_files, a3m_files):
        objects.append(ComPotts_Object.from_hhblits_output(seq_file, a3m_file, **kwargs))
    align_two_objects(objects, output_folder, **kwargs)


def align_one_hot(seq_files, output_folder, **kwargs):
    objects = []
    for seq_file in seq_files :
        objects.append(ComPotts_Object.from_seq_file_to_one_hot(seq_file))
    align_two_objects(objects, output_folder, **kwargs)
