from compotts_wrapper.compute_scores import *
from potts_model import *
import ctypes

COMPOTTS_LOCATION = "./ComPotts"
COMPOTTS_SOLVER = ctypes.CDLL("./compotts_solver.so")
INFINITY = 10000000 


def align_two_potts_models(mrfs, aln_res_file, info_res_file, n_limit_param=INFINITY, iter_limit_param=INFINITY, t_limit=36000, disp_level=1, epsilon=1, v_score_function=scalar_product, w_score_function=scalar_product, gap_open=0, gap_extend=0, w_threshold=0, **kwargs):
    v_scores = compute_v_scores(*mrfs, v_score_function)
    c_v_scores = ctypes.c_void_p(v_scores.ctypes.data)
    edges_maps = [get_edges_map(mrf, w_threshold) for mrf in mrfs]
    w_scores = compute_w_scores(*mrfs, *edges_maps, w_score_function)
    c_w_scores = ctypes.c_void_p(w_scores.ctypes.data)
    c_edges_mapA = ctypes.c_void_p((edges_maps[0].flatten()).ctypes.data)
    c_edges_mapB = ctypes.c_void_p((edges_maps[1].flatten()).ctypes.data)
    selfcompA = compute_selfscore(mrfs[0], v_score_function, w_score_function)
    selfcompB = compute_selfscore(mrfs[1], v_score_function, w_score_function)
    COMPOTTS_SOLVER.call_from_python(c_v_scores, c_w_scores, ctypes.c_int(mrfs[0].ncol), ctypes.c_int(mrfs[1].ncol), c_edges_mapA, c_edges_mapB, ctypes.c_double(selfcompA), ctypes.c_double(selfcompB), ctypes.c_double(gap_open), ctypes.c_double(gap_extend), ctypes.c_char_p(aln_res_file.encode('utf-8')), ctypes.c_char_p(info_res_file.encode('utf-8')), ctypes.c_int(iter_limit_param), ctypes.c_int(t_limit), ctypes.c_int(disp_level), ctypes.c_double(epsilon))


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
