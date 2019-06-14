from compotts_wrapper import compute_scores
from potts_model import *

COMPOTTS_LOCATION = "./ComPotts"


def align_two_potts_models_objects(mrfs, aln_res_file, info_res_file, **kwargs):
    call = COMPOTTS_LOCATION+" --out_aln "+aln_res_file+" --out_info "+info_res_file
    if 'scores_folder' not in kwargs:
        kwargs['scores_folder'] = fm.create_folder('temp_scores/')
    scores_files_dict = compute_scores.get_all_scores_in_files(mrfs, **kwargs)
    for score_option in scores_files_dict:
        call+=" --"+score_option+" "+scores_dict[score_option]
    for key_arg in kwargs:
        if key_arg in possible_compotts_options:
            call+=" --"+key_arg+" "+kwargs[key_arg]
    print("ComPotts call :")
    print(call)
    os.system(call)


def align_two_objects(objects, aln_res_file, info_res_file, **kwargs):
    align_two_potts_models_objects([o.mrf for o in objects], aln_res_file, info_res_file, **kwargs)


def align_two_potts_models(mrf_files, aln_res_file, info_res_file, **kwargs):
    mrf = [Potts_Model.from_msgpack(mrf_file) for mrf_file in mrf_files]
    align_two_potts_models_objects(mrfs, aln_res_file, info_res_file, **kwargs)


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
