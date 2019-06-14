from util import *
import msgpack

def compute_v_scores_dict(mrf1, mrf2, v_score_function, **kwargs):
    v_scores_dict = {}
    for i in range(mrf1.ncol):
        v_scores_dict[i] = {}
        for k in range(mrf2.ncol):
            v_scores_dict[i][k] = v_score_function(mrf1.v[i], mrf2.v[k])
    return v_scores_dict


def compute_self_v_scores(mrf, v_score_function, **kwargs):
    self_v_scores = []
    for i in range(mrf.ncol):
        self_v_scores.append(v_score_function(mrf.v[i], mrf.v[i]))
    return self_v_scores



def compute_self_w_scores(mrf, w_score_function, **kwargs):
    self_w_scores = {}
    for i in range(mrf.ncol):
        self_w_scores[i] = {}
        for j in range(mrf.ncol):
            self_w_scores[i][j] = w_score_function(mrf.w[i][j], mrf.w[i][j])
    return self_w_scores


def compute_w_scores_dict(mrf1, mrf2, w_score_function, **kwargs):
    w_scores_dict = {}
    for i in range(mrf1.ncol-1):
        w_scores_dict[i] = {}
        for j in range(i+1, mrf1.ncol):
            w_scores_dict[i][j] = {}
            for k in range(mrf2.ncol-1):
                w_scores_dict[i][j][k] = {}
                for l in range(k+1, mrf2.ncol):
                    w_scores_dict[i][j][k][l] = w_score_function(mrf1.w[i][j], mrf2.w[k][l])
    return w_scores_dict


def write_dict_in_file(the_dict, filename):
    with open(filename, 'wb') as f:
        f.write(msgpack.dumps(the_dict))
    return filename


def write_list_in_file(the_list, filename):
    with open(filename, 'wb') as f:
        f.write(msgpack.dumps(the_list))
    return filename



def get_all_scores_in_files(mrfs, **kwargs):
    d = {}
    if "v_score_function" in kwargs:
        v_score_function = kwargs["v_score_function"]
    else:
        v_score_function = scalar_product
    if "w_score_function" in kwargs:
        w_score_function = kwargs["w_score_function"]
    else:
        w_score_function = scalar_product
    if 'output_name' in kwargs:
        output_name = kwargs["output_name"]
    else:
        output_name = '_'.join([mrf.name for mrf in mrfs])
    f = kwargs['scores_folder']+output_name
    d["v_scores"]=write_dict_in_file(compute_v_scores_dict(*mrfs, v_score_function, **kwargs), f+"_v_scores")
    d["w_scores"]=write_dict_in_file(compute_w_scores_dict(*mrfs, w_score_function, **kwargs), f+"_w_scores")
    for mrf in mrfs:
        fs = f+"_"+mrf_name
        d[mrf.name+"_self_v_scores"]=write_list_in_file(compute_self_v_scores(mrf, v_score_function, **kwargs), fs+"_self_v_scores")
        d[mrf.name+"_self_w_scores"]=write_dict_in_file(compute_self_w_scores(mrf, w_score_function, **kwargs), fs+"_self_w_scores")
    return d
