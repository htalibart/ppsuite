import matplotlib.pyplot as plt
import seaborn as sns
import inspect

from comutils.global_variables import ALPHABET
from comutils.util import *

from makepotts.potts_model import *

from compotts.compute_scores import *

def get_reordered_v(v, alphabet):
    """ reorders all vi for a given alphabet """
    idx = [ALPHABET.find(a) for a in alphabet]
    return v[:,idx]


def get_reordered_wij(wij, alphabet):
    """ reorders all wij for a given alphabet """
    idx = [ALPHABET.find(a) for a in alphabet]
    new_wij = np.zeros_like(wij)
    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            new_wij[i][j] = wij[idx[i]][idx[j]]
    return new_wij


def plot_heatmap(matrix, center=0, show_figure=True, **kwargs):
    """ plots a heatmap with seaborn """
    plt.figure()
    sns.heatmap(matrix, cmap="RdBu", center=center, **kwargs)
    plt.tick_params(labelsize='xx-small')
    if show_figure:
        plt.show()


def visualize_v_parameters(v, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    tick_space = 3
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]
    v = get_reordered_v(v, alphabet)
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", center=0)
    plt.tick_params(labelsize='xx-small')
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()


def visualize_v_norms(v_norm, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    tick_space = 3
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,len(v_norm))]
    sns.heatmap([v_norm], xticklabels=xticklabels, yticklabels=[], cmap="RdBu", center=0)
    plt.tick_params(labelsize='xx-small')
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()




def visualize_parameters(v, v_norm, w_norm, name, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    """ displays v, ||v|| and ||w|| """
    tick_space = 3
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]

    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,4,1]})
    plt.text(0, 1.5, name)

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

    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()


def visualize_mrf(mrf, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    """ displays MRF parameters """
    visualize_parameters(mrf.v, mrf.get_v_norms(), mrf.get_w_norms(), mrf.name, alphabet=alphabet, start_at_1=start_at_1, show_figure=show_figure)


def visualize_mrf_from_msgpack(msgpack_file, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    mrf = Potts_Model.from_msgpack(msgpack_file)
    visualize_mrf(mrf, alphabet=alphabet, start_at_1=start_at_1, show_figure=show_figure)


def visualize_mrf_difference(mrf1, mrf2, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    v_diff = mrf1.v-mrf2.v
    v_norm_diff = [euclidean_norm(mrf1.v[i])-euclidean_norm(mrf2.v[i]) for i in range(mrf1.ncol)]
    w_norm = mrf1.get_w_norms()-mrf2.get_w_norms()
    name = mrf1.name+"-"+mrf2.name
    visualize_parameters(v_diff, v_norm_diff, w_norm, name, alphabet, start_at_1, show_figure)

    mrf_diff = Potts_Model.from_parameters(mrf1.v-mrf2.v, mrf1.w-mrf2.w, name=name+"_diffmrf")
    visualize_mrf(mrf_diff, alphabet, start_at_1, show_figure)



def visualize_one_sequence(mrf, sequence, show_figure=True):
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
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()


def visualize_v_alignment_from_files(mrf_files, aln_res_file, **kwargs):
    aligned_mrfs = [Potts_Model.from_msgpack(mrf_file) for mrf_file in mrf_files]
    dict_aligned_pos = fm.get_aligned_positions_dict_from_compotts_output_file(aln_res_file)
    visualize_v_alignment(aligned_mrfs, dict_aligned_pos, **kwargs)


def visualize_v_alignment(aligned_mrfs, dict_aligned_pos, alphabet=ALPHABET, start_at_1=True, show_figure=True, tick_space=3):
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1]})
    aligned_pos = list(dict_aligned_pos.values())
    for k in range(2):
        v = get_reordered_v(aligned_mrfs[k].v[aligned_pos[k],:], alphabet)
        xticklabels = [str(aligned_pos[k][i]+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]
        sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[k], center=0)
        ax[k].tick_params(labelsize='xx-small')
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()


def visualize_pos_norms_alignment(aligned_mrfs, dict_aligned_pos, start_at_1=True, show_figure=True, tick_space=3):
    aligned_pos = list(dict_aligned_pos.values())
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1]})
    for k in range(2):
        mrf = aligned_mrfs[k]
        contribs = np.zeros((2, len(aligned_pos[k])))
        for i in range(len(aligned_pos[k])):
            contribs[0][i] = euclidean_norm(mrf.v[aligned_pos[k][i]])
            for j in range(len(aligned_pos[k])):
                contribs[1][i]+= mrf.get_w_norm_at_pos(aligned_pos[k][i],aligned_pos[k][j])
            contribs[1][i]=contribs[1][i]/2
        xticklabels = [str(aligned_pos[k][t]+start_at_1) if (t%tick_space==0) else " " for t in range(0,len(aligned_pos[k]))]
        sns.heatmap(contribs, yticklabels=["v","w"], xticklabels=xticklabels, cmap="RdBu", ax=ax[k], center=0)
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()


def visualize_v_w_scores_at_positions(aligned_mrfs, dict_aligned_pos, show_figure=True, tick_space=3, v_score_function=scalar_product, w_score_function=scalar_product, label_dict=None):
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1]})
    aligned_pos = list(dict_aligned_pos.values())

    if label_dict is None:
        label_dict = dict_aligned_pos
        label_list = aligned_pos
    else:
        label_list = list(label_dict.values())

    xticklabels = ['('+str(label_list[0][k]+start_at_1)+','+str(label_list[1][k]+start_at_1)+')' if (k%tick_space==0) else " " for k in range(len(aligned_pos[0]))]

    v_scores = [v_score_function(aligned_mrfs[0].v[i],aligned_mrfs[1].v[k]) for i,k in zip(aligned_pos[0], aligned_pos[1])]
    sns.heatmap([v_scores], yticklabels=['v'], xticklabels=[], cmap="RdBu", ax=ax[0], center=0)

    w_scores_sums = [sum([w_score_function(aligned_mrfs[0].w[i][j],aligned_mrfs[1].w[k][l]) for j,l in zip(aligned_pos[0], aligned_pos[1])]) for i,k in zip(aligned_pos[0], aligned_pos[1])]
    sns.heatmap([w_scores_sums], yticklabels=['w'], xticklabels=xticklabels, cmap="RdBu", ax=ax[1], center=0)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.draw()
    if show_figure:
        plt.show()



def visualize_v_scores_alignment(aligned_mrfs, dict_aligned_pos, show_figure=True, tick_space=3, label_with_ref=False, label_with_2=False, v_score_function=scalar_product):
    aligned_pos = list(dict_aligned_pos.values())
    plt.figure()
    scores = [v_score_function(aligned_mrfs[0].v[i],aligned_mrfs[1].v[j]) for i,j in zip(aligned_pos[0], aligned_pos[1])]
    if label_with_ref:
        xticklabels = aligned_pos[0]
    elif label_with_2:
        xticklabels = aligned_pos[1]
    else:
        insp = inspect.getargspec(sns.heatmap)
        xticklabels = insp.defaults[insp.args.index('xticklabels')]
    sns.heatmap([scores], yticklabels=["v"], xticklabels=xticklabels, cmap="RdBu", center=0)
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()


def visualize_w_scores_alignment(aligned_mrfs, dict_aligned_pos, show_figure=True, tick_space=3, label_with_ref=False, label_with_2=False, w_score_function=scalar_product, start_at_1=False):
    aligned_pos = list(dict_aligned_pos.values())
    len_aln = len(aligned_pos[0])
    plt.figure()
    scores = np.zeros((len_aln,len_aln))
    for ind_i in range(len_aln):
        for ind_j in range(len_aln):
                scores[ind_i,ind_j] = w_score_function(aligned_mrfs[0].w[aligned_pos[0][ind_i],aligned_pos[0][ind_j]], aligned_mrfs[1].w[aligned_pos[1][ind_i],aligned_pos[1][ind_j]])
    if label_with_ref:
        xticklabels = aligned_pos[0]
        yticklabels = aligned_pos[1]
    elif label_with_2:
        xticklabels = aligned_pos[1]
        yticklabels = aligned_pos[0]
    else:
        #insp = inspect.getargspec(sns.heatmap)
        #xticklabels = insp.defaults[insp.args.index('xticklabels')]
        xticklabels = []
        yticklabels = []
    xticklabels = [str(xi+start_at_1) if (xi%tick_space==0) else " " for xi in xticklabels]
    yticklabels = [str(xi+start_at_1) if (xi%tick_space==0) else " " for xi in yticklabels]
    sns.heatmap(scores, yticklabels=yticklabels, xticklabels=xticklabels, cmap="RdBu", center=0)
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()



def visualize_v_w_scores_alignment(aligned_mrfs, dict_aligned_pos, show_figure=True, tick_space=3, v_score_function=scalar_product, w_score_function=scalar_product, alphabet=ALPHABET, start_at_1=False, label_dict=None):
    aligned_pos = list(dict_aligned_pos.values())

    if label_dict is None:
        label_dict = dict_aligned_pos
        label_list = aligned_pos
    else:
        label_list = list(label_dict.values())


    len_aln = len(aligned_pos[0])

    fig, ax = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1,6,1]})

    # v alignment : sqrt(vi(a)*vk(a))
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
    v_scores = [v_score_function(aligned_mrfs[0].v[i],aligned_mrfs[1].v[j]) for i,j in zip(aligned_pos[0], aligned_pos[1])]
    sns.heatmap([v_scores], xticklabels=[], yticklabels=['v'], cmap="RdBu", center=0, ax=ax[1])
    #ax[2].tick_params(labelsize='xx-small')

    # w scores
    w_scores = np.zeros((len_aln,len_aln))
    for ind_i in range(len_aln):
        for ind_j in range(len_aln):
                w_scores[ind_i,ind_j] = w_score_function(aligned_mrfs[0].w[aligned_pos[0][ind_i],aligned_pos[0][ind_j]], aligned_mrfs[1].w[aligned_pos[1][ind_i],aligned_pos[1][ind_j]])
    xticklabels = [(label_list[0][k]+start_at_1,label_list[1][k]+start_at_1) for k in range(len(aligned_pos[0]))]
    xticklabels = [xi if (i%tick_space==0) else " " for i, xi in enumerate(xticklabels)]
    yticklabels=xticklabels
    sns.heatmap(w_scores, xticklabels=[], yticklabels=yticklabels, cmap="RdBu", center=0, ax=ax[2])
    ax[2].tick_params(labelsize='x-small')


    # w scores contributions
    w_scores_sums = [sum([w_score_function(aligned_mrfs[0].w[i][j],aligned_mrfs[1].w[k][l]) for j,l in zip(aligned_pos[0], aligned_pos[1])]) for i,k in zip(aligned_pos[0], aligned_pos[1])]
    sns.heatmap([w_scores_sums], yticklabels=['w'], xticklabels=xticklabels, cmap="RdBu", ax=ax[3], center=0)
    ax[3].tick_params(labelsize='x-small')


    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.draw()
    if show_figure:
        plt.show()
    return fig



def visualize_v_norm_alignment(aligned_mrfs, dict_aligned_pos, start_at_1=True, show_figure=True, tick_space=3):

    aligned_pos = list(dict_aligned_pos.values())
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,1]})
    for k in range(2):
        v_norms = [euclidean_norm(aligned_mrfs[k].v[i]) for i in aligned_pos[k]]
        xticklabels = [str(aligned_pos[k][i]+start_at_1) if (i%tick_space==0) else " " for i in range(len(v_norms))]
        sns.heatmap([v_norms], xticklabels=xticklabels, cmap="RdBu", ax=ax[k], center=0)
        ax[k].tick_params(labelsize='xx-small')
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()



def visualize_letters_alignment(aligned_mrfs, aligned_pos, alphabet=ALPHABET, start_at_1=True, show_figure=True, tick_space=3):
    plt.figure()

    vs = [get_reordered_v(aligned_mrfs[k].v[aligned_pos[k],:], alphabet) for k in range(2)]

    v_align = np.zeros((len(alphabet)*2, len(vs[0])))
    for a in range(len(alphabet)):
        for k in range(2):
            for i in range(len(vs[0])):
                v_align[2*a+k][i] = vs[k][i][a]

    sns.heatmap(v_align, yticklabels=[alphabet[j//2] if (j%2==0) else " " for j in range(2*len(alphabet))], xticklabels=[], cmap="RdBu", center=0)
    plt.tick_params(labelsize='xx-small')
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()



def plot_one_vi(vi, alphabet=ALPHABET, **kwargs):
    if alphabet in kwargs:
        idx = [ALPHABET.find(a) for a in alphabet]
        new_vi = vi[idx]
    else:
        new_vi = vi
    plt.figure()
    sns.heatmap(vi.reshape(vi.shape[0], 1), square=True, cmap="RdBu", center=0, yticklabels=alphabet, xticklabels=[], **kwargs)
    plt.margins(0,0)
    plt.tight_layout()
    plt.show()


def plot_one_wij(wij, alphabet=ALPHABET, center=0, show_figure=True, **kwargs):
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
