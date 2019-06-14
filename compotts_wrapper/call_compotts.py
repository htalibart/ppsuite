from compotts_wrapper.compute_scores import *
from potts_model import *
import ctypes

COMPOTTS_LOCATION = "./ComPotts"
COMPOTTS_SOLVER = ctypes.CDLL("./compotts_solver.so")
INFINITY = 10000000 


def align_two_potts_models(mrfs, aln_res_file, info_res_file, n_limit_param=INFINITY, iter_limit_param=INFINITY, t_limit=36000, disp_level=1, epsilon=1, v_score_function=scalar_product, w_score_function=scalar_product, **kwargs):
    v_scores = compute_v_scores(*mrfs, v_score_function)
    c_v_scores = ctypes.c_void_p(v_scores.ctypes.data)
    w_scores = compute_w_scores(*mrfs, w_score_function)
    c_w_scores = ctypes.c_void_p(w_scores.ctypes.data)
    edges_mapA = get_edges_map(mrfs[0])
    c_edges_mapA = ctypes.c_void_p(edges_mapA.ctypes.data)
    edges_mapB = get_edges_map(mrfs[1])
    c_edges_mapB = ctypes.c_void_p(edges_mapB.ctypes.data)
    selfcompA = compute_selfscore(mrfs[0])
    selfcompB = compute_selfscore(mrfs[1])
    COMPOTTS_SOLVER.call_from_python(c_v_scores, c_w_scores, ctypes.c_int(mrfs[0].ncol), ctypes.c_int(mrfs[1].ncol), c_edges_mapA, c_edges_mapB, ctypes.c_double(selfcompA), ctypes.c_double(selfcompB), ctypes.c_double(gap_open), ctypes.c_double(gap_extend), ctypes.c_char_p(aln_res_file), ctypes.c_char_p(info_res_file), ctypes.c_int(iter_limit_param), ctypes.c_int(t_limit), ctypes.c_int(disp_level), ctypes.c_double(c_epsilon))


def align_two_objects(objects, aln_res_file, info_res_file, **kwargs):
    align_two_potts_models([o.mrf for o in objects], aln_res_file, info_res_file, **kwargs)


def align_hhblits_output(seq_files, a3m_files, aln_res_file, info_res_file, **kwargs):
    objects = []
    for seq_file, a3m_file in zip(seq_files, a3m_files):
        objects.append(compotts_object.from_hhblits_output(seq_file, a3m_file, kwargs))
    align_two_objects(objects, aln_res_file, info_res_file, kwargs)


def align_one_hot(seq_files, aln_res_file, info_res_file, **kwargs):
    objects = []
    for seq_file in seq_files :
        objects.append(compotts_object.from_seq_file_to_one_hot(seq_file))


# def multiple_alignment() TODO
