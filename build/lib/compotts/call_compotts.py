import ctypes

from compotts.compute_scores import *
from compotts.compotts_object import *
from basic_modules.potts_model import *
import basic_modules.files_management as fm

COMPOTTS_SOLVER = fm.get_compots_solver()
INFINITY = 10000000 

# TODO stocker les scores en utilisant la symétrie pour prendre moins de mémoire
# TODO aln_res_file et info_res_file comme sortie de ComPotts -> variables

def align_two_potts_models(mrfs, output_folder, n_limit_param=INFINITY, iter_limit_param=INFINITY, t_limit=36000, disp_level=1, epsilon=1, v_score_function=scalar_product, w_score_function=scalar_product, gap_open=0, gap_extend=0, w_threshold=0, use_w=True, **kwargs):

    aln_res_file = fm.get_aln_res_file_name(output_folder)
    info_res_file = fm.get_info_res_file_name(output_folder)

    v_scores = np.ascontiguousarray(compute_v_scores(*mrfs, v_score_function).flatten())
    c_v_scores = ctypes.c_void_p(v_scores.ctypes.data)

    if use_w:
        edges_maps = [get_edges_map(mrf, w_threshold) for mrf in mrfs]
        w_scores = compute_w_scores(*mrfs, *edges_maps, w_score_function)
    else:
        edges_maps = [np.zeros((mrf.w.shape[0:2])) for mrf in mrfs]
        w_scores = np.zeros((mrfs[0].w.shape[0],mrfs[0].w.shape[1],mrfs[1].w.shape[0],mrfs[1].w.shape[1]))
    
    c_w_scores = ctypes.c_void_p(w_scores.ctypes.data)

    c_int_p = ctypes.POINTER(ctypes.c_int)
    c_edges_maps = [np.ascontiguousarray(edges_map.flatten(), dtype=np.int32).ctypes.data_as(c_int_p) for edges_map in edges_maps]

    selfcomps = [compute_selfscore(mrf, edges_map, v_score_function, w_score_function) for mrf, edges_map in zip(mrfs, edges_maps)]
    
    COMPOTTS_SOLVER.call_from_python(c_v_scores, c_w_scores, *[ctypes.c_int(mrf.ncol) for mrf in mrfs], *c_edges_maps, *[ctypes.c_double(selfcomp) for selfcomp in selfcomps], ctypes.c_double(gap_open), ctypes.c_double(gap_extend), ctypes.c_char_p(aln_res_file.encode('utf-8')), ctypes.c_char_p(info_res_file.encode('utf-8')), ctypes.c_int(iter_limit_param), ctypes.c_int(t_limit), ctypes.c_int(disp_level), ctypes.c_double(epsilon))

    aligned_positions_dict = fm.get_aligned_positions_dict_from_compotts_output_file(aln_res_file)
    infos_solver = fm.get_infos_solver_dict_from_compotts_output_file(info_res_file)

    return aligned_positions_dict, infos_solver



def align_two_potts_models_from_files(mrf_files, output_folder, **kwargs):
    mrfs = [Potts_Model.from_msgpack(mrf_file) for mrf_file in mrf_files]
    return align_two_potts_models(mrfs, output_folder, **kwargs)


def align_two_objects(objects, output_folder, **kwargs):
    return align_two_potts_models([o.mrf for o in objects], output_folder, **kwargs)

