import ctypes
import time
import pandas as pd

from compotts.compute_scores import *
from makepotts.potts_object import *
from makepotts.potts_model import *
import comutils.files_management as fm

import pkg_resources
COMPOTTS_CPP_LIBRARY = pkg_resources.resource_filename('compotts', 'compotts_solver.so')

COMPOTTS_SOLVER = ctypes.CDLL(COMPOTTS_CPP_LIBRARY)
INFINITY = 1000000000

def align_two_potts_models(mrfs, output_folder, n_limit_param=INFINITY, iter_limit_param=1000, t_limit=36000, disp_level=1, precision_method="arbitrary_1", w_threshold_method="none", use_w=True, use_v=True, gamma=1.0, theta=0.9, stepsize_min=0.000000005, nb_non_increasing_steps_max=500, alpha_w=1, gap_open=8, gap_extend=0, sim_min=0.1, **kwargs):

    aln_res_file = fm.get_aln_res_file_name(output_folder)
    info_res_file = fm.get_info_res_file_name(output_folder)

    time_start = time.time()
    
    
    c_float_p = ctypes.POINTER(ctypes.c_float) # defining type for ctypes : pointer to float

    # handle v
    if use_v:
        v_flats = [np.ascontiguousarray(mrf.v.flatten()) for mrf in mrfs]
    else:
        v_flats = [np.ascontiguousarray(np.zeros(mrf.v.shape)).flatten() for mrf in mrfs]
    c_vs = [vflat.astype(np.float32).ctypes.data_as(c_float_p) for vflat in v_flats]

    if use_w:
        edges_maps = [get_edges_map(mrf, w_threshold_method) for mrf in mrfs]
        w_flats = [np.ascontiguousarray(mrf.w.flatten()) for mrf in mrfs]
    else:
        edges_maps = [np.zeros((mrf.w.shape[0:2])) for mrf in mrfs]
        w_flats = [np.ascontiguousarray(np.zeros(mrf.w.shape)).flatten() for mrf in mrfs]
    c_ws = [wflat.astype(np.float32).ctypes.data_as(c_float_p) for wflat in w_flats]
        
    selfcomps = [compute_selfscore(mrf, edges_map, use_v, use_w, alpha_w, **kwargs) for mrf, edges_map in zip(mrfs, edges_maps)]
    epsilon = get_epsilon(precision_method, selfcomps)

    print("gap open=",gap_open)

    c_int_p = ctypes.POINTER(ctypes.c_int)
    c_edges_maps = [np.ascontiguousarray(edges_map.flatten(), dtype=np.int32).ctypes.data_as(c_int_p) for edges_map in edges_maps]

    score_min = (1/2)*sim_min*sum(selfcomps); 


    COMPOTTS_SOLVER.call_from_python.argtypes=[c_float_p, c_float_p, c_float_p, c_float_p, ctypes.c_int, ctypes.c_int, c_int_p, c_int_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double]

    COMPOTTS_SOLVER.call_from_python(*c_vs, *c_ws, *[ctypes.c_int(mrf.ncol) for mrf in mrfs], *c_edges_maps, *[ctypes.c_double(selfcomp) for selfcomp in selfcomps], ctypes.c_double(gap_open), ctypes.c_double(gap_extend), ctypes.c_char_p(str(aln_res_file).encode('utf-8')), ctypes.c_char_p(str(info_res_file).encode('utf-8')), ctypes.c_int(n_limit_param), ctypes.c_int(iter_limit_param), ctypes.c_double(t_limit), ctypes.c_int(disp_level), ctypes.c_double(epsilon), ctypes.c_double(gamma), ctypes.c_double(theta), ctypes.c_double(stepsize_min), ctypes.c_int(nb_non_increasing_steps_max), ctypes.c_double(score_min))

    total_computation_time = time.time()-time_start

    df = pd.read_csv(info_res_file)
    df['total_compotts_time'] = total_computation_time
    df.to_csv(info_res_file, index=False, na_rep='nan')
    infos_solver = fm.get_infos_solver_dict_from_compotts_output_file(info_res_file)
    
    if not math.isnan(df['similarity_global']):
        aligned_positions_dict = fm.get_aligned_positions_dict_from_compotts_output_file(aln_res_file)
    else:
        aligned_positions_dict = {}

    return aligned_positions_dict, infos_solver


def align_two_objects(objects, output_folder, **kwargs):
    return align_two_potts_models([o.potts_model for o in objects], output_folder, **kwargs)

