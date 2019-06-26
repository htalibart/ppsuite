import matplotlib.pyplot as plt
import seaborn as sns

from basic_modules.global_variables import ALPHABET
from basic_modules.util import *

def get_reordered_v(v, alphabet):
    idx = [ALPHABET.find(a) for a in alphabet]
    return v[:,idx]


def get_reordered_wij(wij, alphabet):
    idx = [ALPHABET.find(a) for a in alphabet]
    new_wij = np.zeros_like(wij)
    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            new_wij[i][j] = wij[idx[i]][idx[j]]
    return new_wij


def plot_heatmap(v_scores, center=0, show_figure=True, **kwargs):
    plt.figure()
    sns.heatmap(np.transpose(v_scores), cmap="RdBu", center=center, **kwargs)
    plt.tick_params(labelsize='xx-small')
    if show_figure:
        plt.show()


def visualize_parameters(v, v_norm, w_norm, name, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    """ affiche les paramètres v, ||v|| et ||w|| """
    tick_space = 3
    xticklabels = [str(i+start_at_1) if (i%tick_space==0) else " " for i in range(0,v.shape[0])]

    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False, gridspec_kw={'height_ratios':[1,4,1]})
    plt.text(0, 1.5, name)

    v = get_reordered_v(v, alphabet)
    sns.heatmap(np.transpose(v), yticklabels=alphabet, xticklabels=xticklabels, cmap="RdBu", ax=ax[0], center=0)
    ax[0].tick_params(labelsize='xx-small')

    sns.heatmap(w_norm, xticklabels=xticklabels, yticklabels=xticklabels, cmap="RdBu", center=0, ax=ax[1])
    ax[1].tick_params(labelsize='xx-small')

    sns.heatmap([v_norm], xticklabels=xticklabels, yticklabels=[], cmap="RdBu", center=0, ax=ax[2])
    ax[2].tick_params(labelsize='xx-small')

    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()



def visualize_mrf(mrf, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    visualize_parameters(mrf.v, mrf.get_v_norms(), mrf.get_w_norms(), mrf.name, alphabet, start_at_1, show_figure)



def visualize_mrf_difference(mrf1, mrf2, alphabet=ALPHABET, start_at_1=True, show_figure=True):
    """ visualiser la différence entre deux MRFs """
    v_diff = mrf1.v-mrf2.v
    v_norm_diff = [euclidean_norm(mrf1.v[i])-euclidean_norm(mrf2.v[i]) for i in range(mrf1.ncol)]
    w_norm = mrf1.get_w_norms()-mrf2.get_w_norms()
    name = mrf1.name+"-"+mrf2.name
    visualize_parameters(v_diff, v_norm_diff, w_norm, alphabet, name, start_at_1, show_figure)

    mrf_diff = Potts_Model.from_parameters(mrf1.ncol, mrf1.v-mrf2.v, mrf1.w-mrf2.w, name+"_diffmrf")
    visualize_mrf(mrf_diff, alphabet, start_at_1, show_figure)



def visualize_one_sequence(mrf, sequence, show_figure=True):
    """ afficher les paramètres v_i_a et w_ij_ab du mrf pour les a et b correspondant à la séquence """
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



def visualize_pos_contributions_alignment(aligned_mrfs, aligned_pos, start_at_1=True, show_figure=True, tick_space=3):
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


def visualize_v_w_scores_alignment(aligned_mrfs, aligned_pos, aligned_v_scores, aligned_w_scores, show_figure=True, tick_space=3):
    plt.figure()
    scores = np.zeros((2, len(aligned_pos[0])))
    scores[0] = aligned_v_scores
    scores[1] = aligned_w_scores
    sns.heatmap(scores, yticklabels=["v","w"], cmap="RdBu", center=0)
    plt.tight_layout()
    plt.draw()
    if show_figure:
        plt.show()



def visualize_v_norm_alignment(aligned_mrfs, aligned_pos, start_at_1=True, show_figure=True, tick_space=3):
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
    plt.show()
