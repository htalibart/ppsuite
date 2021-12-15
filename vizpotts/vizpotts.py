import matplotlib.pyplot as plt
import seaborn as sns
import inspect

from comutils.global_variables import ALPHABET
from comutils.util import *

from makepotts.potts_model import *

from ppalign.compute_scores import *

def end_visual(tight_layout=True, show_figure=True, **kwargs):
    """ factorizing matplotlib options always used """
    if tight_layout:
        plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()

def get_reordered_v(v, alphabet):
    """ reorders all vi for a given alphabet """
    q = v.shape[1]
    idx = [ALPHABET[:q].find(a) for a in alphabet[:q]]
    return v[:,idx]


def get_reordered_wij(wij, alphabet):
    """ reorders all wij for a given alphabet """
    q = wij.shape[0]
    idx = [ALPHABET[:q].find(a) for a in alphabet[:q]]
    new_wij = np.zeros_like(wij)
    for a in range(len(alphabet)):
        for b in range(len(alphabet)):
            new_wij[a][b] = wij[idx[a]][idx[b]]
    return new_wij


def plot_heatmap(matrix, center=0, **kwargs):
    """ plots a heatmap with seaborn """
    plt.figure()
    sns.heatmap(matrix, cmap="RdBu", center=center, **kwargs)
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)


def visualize_insertion_penalties(insertion_penalties, **kwargs):
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1]})
    sns.heatmap([insertion_penalties["open"]], ax=ax[0])
    ax[0].collections[0].colorbar.set_label("gap open")
    sns.heatmap([insertion_penalties["extend"]], ax=ax[1])
    ax[1].collections[0].colorbar.set_label("gap extend")
    end_visual(**kwargs)


def visualize_v_parameters(v, alphabet=ALPHABET, start_at_1=True, tick_space=3, figsize=(10,2), center=0, **kwargs):
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]
    v = get_reordered_v(v, alphabet)
    plt.figure(figsize=figsize)
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", center=center, cbar_kws={'label': r'$v_i(a)$'})
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)


def visualize_v_norms(v_norm, start_at_1=True, tick_space=3, figsize=(10,2), colorbar_label=r'$||v_i||$', **kwargs):
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,len(v_norm))]
    plt.figure(figsize=figsize)
    sns.heatmap([v_norm], xticklabels=xticklabels, yticklabels=[], cmap="RdBu", center=0, cbar_kws={'label': colorbar_label})
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)
  

def visualize_w_norms(w_norm, start_at_1=True, tick_space=3, figsize=(10,9), tight_layout=True, colorbar_label = r'$||w_{ij}||$', **kwargs):
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,len(w_norm))]
    plt.figure(figsize=figsize)
    sns.heatmap(w_norm, xticklabels=xticklabels, yticklabels=xticklabels, cmap="RdBu", center=0, cbar_kws={'label': colorbar_label})
    plt.tick_params(labelsize='xx-small')
    end_visual(**kwargs)



def visualize_parameters(v, v_norm, w_norm, name, alphabet=ALPHABET, start_at_1=True, **kwargs):
    """ displays v, ||v|| and ||w|| """
    tick_space = 3
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]

    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,4,1]})
    plt.text(0, 2.5, '...'+name[-50:])

    v = get_reordered_v(v, alphabet)
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[0], center=0)
    ax[0].tick_params(labelsize='xx-small')
    ax[0].set_xlabel('i')
    ax[0].set_ylabel('a')
    ax[0].collections[0].colorbar.set_label("vi(a)")
    #ax[0].text(-5,10,'v',fontsize=15)

    sns.heatmap(w_norm, xticklabels=xticklabels, yticklabels=xticklabels, cmap="RdBu", center=0, ax=ax[1])
    ax[1].tick_params(labelsize='xx-small')
    ax[1].set_xlabel('i')
    ax[1].set_ylabel('j')
    ax[1].collections[0].colorbar.set_label("||wij||")

    sns.heatmap([v_norm], xticklabels=xticklabels, yticklabels=[], cmap="RdBu", center=0, ax=ax[2])
    ax[2].tick_params(labelsize='xx-small')
    ax[2].set_xlabel('i')
    ax[2].collections[0].colorbar.set_label("||vi||")

    end_visual(**kwargs)



def visualize_mrf(mrf, alphabet=ALPHABET, start_at_1=True, **kwargs):
    """ displays MRF parameters """
    visualize_parameters(mrf.v, mrf.get_v_norms(), mrf.get_w_norms(), mrf.name, alphabet=alphabet, start_at_1=start_at_1, **kwargs)



def visualize_one_sequence(mrf, sequence, **kwargs):
    """ visualization of parameters v_i_a and w_ij_ab of Potts model @mrf for a and b in sequence @sequence """
    xticklabels = [s for s in sequence]
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,6]})
    v_seq = [mrf.v[i,code(sequence[i])] for i in range(len(sequence))]
    sns.heatmap([v_seq], cmap="RdBu", ax=ax[0], xticklabels=xticklabels, center=0)

    w_seq = np.zeros((len(sequence), len(sequence)))
    for i in range(len(sequence)):
        for j in range(len(sequence)):
            w_seq[i,j] = mrf.w[i,j,mrf.code(sequence[i]),mrf.code(sequence[j])]
    sns.heatmap(w_seq, cmap="RdBu", ax=ax[1], xticklabels=xticklabels, yticklabels=xticklabels, center=0)
    end_visual(**kwargs)



def visualize_v_alignment(aligned_mrfs, aln_res_file, alphabet=ALPHABET, start_at_1=True, tick_space=3, label_dict=None, **kwargs):
    dict_aligned_pos = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_res_file)
    aligned_pos = [dict_aligned_pos["pos_ref"], dict_aligned_pos["pos_2"]]
    if label_dict is None:
        label_dict = dict_aligned_pos
        label_list = aligned_pos
    else:
        label_list = list(label_dict.values())

    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1]})
    for k in range(2):
        v = get_reordered_v(aligned_mrfs[k].v[aligned_pos[k],:], alphabet)
        xticklabels = [str(label_list[k][i]+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]
        sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[k], center=0)
        ax[k].tick_params(labelsize='xx-small')
    end_visual(**kwargs)
   

def visualize_v_alignment_with_scalar_product(aligned_mrfs, aln_res_file, alphabet=ALPHABET, start_at_1=True, tick_space=3, label_dict=None, v_score_function=scalar_product, **kwargs):
    dict_aligned_pos = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_res_file)
    aligned_pos = [dict_aligned_pos["pos_ref"], dict_aligned_pos["pos_2"]]
    if label_dict is None:
        label_dict = dict_aligned_pos
        label_list = aligned_pos
    else:
        label_list = list(label_dict.values())

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1,1,1]})

    # vi alignment
    for k in range(2):
        v = get_reordered_v(aligned_mrfs[k].v[aligned_pos[k],:], alphabet)
        xticklabels = [str(label_list[k][i]+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]
        sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[k], center=0)
        ax[k].tick_params(labelsize='xx-small')

    # vi*vk
    len_aln = len(aligned_pos[0])
    len_v = len(aligned_mrfs[0].v[0])
    letter_v_scores = np.zeros((len_aln, len_v))
    for ind_i in range(len_aln):
        for a in range(len_v):
            letter_v_scores[ind_i][a] = aligned_mrfs[0].v[aligned_pos[0][ind_i]][a]*aligned_mrfs[1].v[aligned_pos[1][ind_i]][a]
    v = get_reordered_v(letter_v_scores, alphabet)
    xticklabels = []
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[2], center=0)
    ax[2].tick_params(labelsize='xx-small')

    # v scores
    v_scores = get_v_scores_for_alignment(aligned_mrfs, dict_aligned_pos, v_score_function=v_score_function, **kwargs)
    sns.heatmap([v_scores], xticklabels=[], yticklabels=['v'], cmap="RdBu", center=0, ax=ax[3])


    end_visual(**kwargs)



def visualize_v_w_scores_alignment(aligned_mrfs, aln_res_file, show_figure=True, tick_space=3, v_score_function=scalar_product, w_score_function=scalar_product, alpha_w=1, alphabet=ALPHABET, start_at_1=False, label_dict=None, print_text=True, **kwargs):
    dict_aligned_pos = fm.get_aligned_positions_dict_from_ppalign_output_file(aln_res_file)
    aligned_pos = [dict_aligned_pos["pos_ref"], dict_aligned_pos["pos_2"]]

    if label_dict is None:
        label_dict = dict_aligned_pos
        label_list = aligned_pos
    else:
        label_list = list(label_dict.values())


    len_aln = len(aligned_pos[0])

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1,1,6]})

    # v alignment : vi(a)*vk(a)
    len_v = len(aligned_mrfs[0].v[0])
    letter_v_scores = np.zeros((len_aln, len_v))
    for ind_i in range(len_aln):
        for a in range(len_v):
            letter_v_scores[ind_i][a] = aligned_mrfs[0].v[aligned_pos[0][ind_i]][a]*aligned_mrfs[1].v[aligned_pos[1][ind_i]][a]
    v = get_reordered_v(letter_v_scores, alphabet)
    xticklabels = []
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[0], center=0)
    ax[0].tick_params(labelsize='xx-small')


    # v scores alignment
    v_scores = get_v_scores_for_alignment(aligned_mrfs, dict_aligned_pos, v_score_function=v_score_function, **kwargs)
    sns.heatmap([v_scores], xticklabels=[], yticklabels=['v'], cmap="RdBu", center=0, ax=ax[1])
    #ax[2].tick_params(labelsize='xx-small')

    # w scores contributions
    w_scores_sums = [sum([get_wij_wkl_score(aligned_mrfs[0].w[i][j],aligned_mrfs[1].w[k][l], w_score_function=w_score_function, **kwargs) for j,l in zip(aligned_pos[0], aligned_pos[1])]) for i,k in zip(aligned_pos[0], aligned_pos[1])]
    sns.heatmap([w_scores_sums], yticklabels=['w'], xticklabels=[], cmap="RdBu", ax=ax[2], center=0)
    ax[2].tick_params(labelsize='x-small')


    # w scores
    w_scores = get_w_scores_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs)
    xticklabels = [(label_list[0][k]+start_at_1,label_list[1][k]+start_at_1) for k in range(len(aligned_pos[0]))]
    xticklabels = [xi if (i%tick_space==0) else " " for i, xi in enumerate(xticklabels)]
    yticklabels=xticklabels
    sns.heatmap(w_scores, xticklabels=xticklabels, yticklabels=yticklabels, cmap="RdBu", center=0, ax=ax[3])
    ax[3].tick_params(labelsize='x-small')



    # print scores
    if print_text:
        text="total PPalign score : "+"{:4.4f}".format(get_score_for_alignment(aligned_mrfs, dict_aligned_pos, alpha_w=alpha_w, **kwargs))+"\npositional score : "+"{:4.4f}".format(get_v_score_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs))+"\ncoupling score : "+"{:4.4f}".format(get_w_score_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs))+" x "+str(alpha_w)+"="+"{:4.4f}".format(get_w_score_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs)*alpha_w)+"\ngap cost : "+str(get_total_gap_cost(dict_aligned_pos, sequence_lengths=[mrf.ncol for mrf in aligned_mrfs], gap_open=kwargs["gap_open"]))+"\nnb positions aligned : "+str(len_aln)
        #text="positional score : "+"{:4.4f}".format(get_v_score_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs))+"\ncoupling score : "+"{:4.4f}".format(get_w_score_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs))+" x "+str(alpha_w)+"="+"{:4.4f}".format(get_w_score_for_alignment(aligned_mrfs, dict_aligned_pos, **kwargs)*alpha_w)+"\nnb positions aligned : "+str(len_aln)
        plt.gca()
        plt.subplots_adjust(bottom=0.5)
        plt.figtext(0.05,0.01, text, fontsize=10, va="bottom", ha="left")

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.draw()
    if show_figure:
        plt.show()
    return fig



def plot_one_vi(vi, alphabet=ALPHABET, show_figure=True, **kwargs):
    if alphabet in kwargs:
        idx = [ALPHABET.find(a) for a in alphabet]
        new_vi = vi[idx]
    else:
        new_vi = vi
    plt.figure()
    sns.heatmap(vi.reshape(vi.shape[0], 1), square=True, cmap="RdBu", center=0, yticklabels=alphabet, xticklabels=[], **kwargs)
    plt.margins(0,0)
    plt.tight_layout()
    if show_figure:
        plt.show()


def plot_one_wij(wij, alphabet=ALPHABET, center=0, show_figure=True, **kwargs):
    q = wij.shape[0]
    alphabet=alphabet[:q]
    wij = get_reordered_wij(wij, alphabet)
    plt.figure()
    sns.heatmap(wij, cmap="RdBu", center=center, xticklabels=alphabet, yticklabels=alphabet, **kwargs)
    plt.tick_params(labelsize='xx-small')
    if show_figure:
        plt.show()


def plot_v_scores(mrfs, v_score_function=scalar_product, **kwargs):
    v_scores = compute_v_scores(mrfs[0], mrfs[1], v_score_function, **kwargs)
    plot_heatmap(v_scores, **kwargs)


def plot_v_scores_path(mrfs, dict_aligned_positions, v_score_function=scalar_product, **kwargs):
    v_scores = compute_v_scores(mrfs[0], mrfs[1], v_score_function, **kwargs)
    v_scores_path = np.zeros_like(v_scores)
    for pos in range(len(dict_aligned_positions['pos_ref'])):
        i=dict_aligned_positions['pos_ref'][pos]
        k=dict_aligned_positions['pos_2'][pos]
        v_scores_path[i][k] = v_scores[i][k]
    plot_heatmap(v_scores_path, **kwargs)
    print("sum =",np.sum(v_scores_path))
