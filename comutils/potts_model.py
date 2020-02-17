import os
import pathlib
import numpy as np
import numpy.linalg as LA
import math
import msgpack
import json
import subprocess

from comutils.util import *
from comutils import files_management as fm
from comutils import pseudocounts

POSSIBLE_CCMPRED_OPTIONS = ["wt-simple", "wt-simple", "wt-uniform", "wt-cutoff", "reg-lambda-single", "reg-lambda-pair-factor", "reg-L2", "reg-noscaling", "reg-scale-by-L", "v-center", "v-zero", "max-gap-pos", "max-gap_seq", "pc-uniform", "pc-submat", "pc-constant", "pc-none", "pc-single-count", "pc-pair-count", "maxit", "ofn-pll", "ofn-cd", "pc-pair-submat", "persistent", "no-decay", "nr-markov-chains"] # TODO mettre toutes les options

class Potts_Model:

    def __init__(self, v, w, **kwargs):
        self.v = v
        self.ncol = len(self.v)
        self.w = w
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = "Billy"
        print("w norm: ", self.get_w_norm(), "normalized w norm: ", self.get_normalized_w_norm())


    @classmethod
    def from_msgpack(cls, binary_file, **kwargs):
        """
            initialize MRF from msgpack file given by CCMpredPy
        """
        with open(str(binary_file), 'rb') as data_file:
            df = msgpack.unpackb(data_file.read())
            """
            df est un dictionnaire :
                b'format'
                b'ncol'
                b'x_single': np.array(ncol, 20)) 
                b'x_pair' : np.array(ncol, ncol, 21, 21)
                b'meta'
            """
        print("getting Potts model from "+str(binary_file))
        fm.check_if_file_ok(binary_file)
        ncol = df[b'ncol']
        v_20 = np.array(df[b'x_single']).reshape((ncol,20))
        v = np.zeros((ncol,21))
        v[:,:-1] = v_20
        w = np.zeros((ncol, ncol, 21, 21))
        for p in df[b'x_pair'].values():
            i = p[b'i']
            j = p[b'j']
            mat = np.array(p[b'x']).reshape((21, 21))
            w[i, j, :, :] = mat
            w[j, i, :, :] = mat.T
        if 'name' not in kwargs:
            kwargs['name'] = str(binary_file).replace('/','-')
            if kwargs['name'].startswith('-'):
                kwargs['name'] = kwargs['name'][1:]
        mrf = cls(v, w, **kwargs)
        mrf.binary_file = binary_file
        return mrf


    @classmethod
    def from_training_set(cls, aln_file, binary_file, write_readme=True, readme_file=None, **kwargs):
        """
            initialize MRF from training set
        """
        fm.check_if_file_ok(aln_file)
        call = "ccmpred "+str(aln_file)+ " -b "+str(binary_file)
        for key_arg in kwargs:
            arg_ccm = key_arg.replace('_', '-')
            if arg_ccm in POSSIBLE_CCMPRED_OPTIONS:
                if isinstance(kwargs[key_arg], bool):
                    if kwargs[key_arg] is True:
                        arg_value=""
                        call+=" --"+arg_ccm+" "+arg_value
                else:
                    arg_value=str(kwargs[key_arg])
                    call+=" --"+arg_ccm+" "+arg_value
        if write_readme:
            if readme_file is None:
                readme_file = pathlib.Path(str(binary_file)[:-len(".mrf")]+"_mrf_README.txt")
            with readme_file.open(mode='w') as f:
                json.dump(call, f, default=str)
        print(call)
        subprocess.Popen(call, shell=True).wait()
        if not os.path.exists(str(binary_file)):
            raise Exception("CCMpredPy wasn't able to infer the MRF. Protein is probably too long ?")
        mrf = cls.from_msgpack(binary_file)
        mrf.training_set = pathlib.Path(aln_file)
        return mrf


    @classmethod
    def from_sequence_file_to_one_hot(cls, seq_file, **kwargs):
        """ one hot encoding """
        fm.check_if_file_ok(seq_file)
        seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        x = code_whole_seq(seq)
        v = np.zeros((len(x),q))
        for i in range(len(x)):
            v[i,x[i]]=1
        w = np.zeros((len(x),len(x),q,q,)) # TODO laisser à 0 ? # , ?
        for i in range(len(x)):
            for j in range(len(x)):
                w[i,j,x[i],x[j]] = 1
        obj = cls.from_parameters(v, w, **kwargs)
        obj.training_set = seq_file
        return obj



    @classmethod
    def from_sequence_file_with_submat(cls, seq_file, npc=1, **kwargs):
        """ substitution matrix pseudocounts """
        fm.check_if_file_ok(seq_file)
        seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        tau = npc/(1+npc)
        x = code_whole_seq(seq)
        v = np.zeros((len(x), q))
        v = np.zeros((len(x),q))

        tau = npc/(1+npc)
        for i in range(len(x)):
            fi = np.zeros(q-1)
            log_sum=0
            for a in range(q-1):
                fi[a] = (1-tau)*(a==x[i]) + tau*pseudocounts.get_cond_proba(a,x[i])
                log_sum+=math.log(fi[a])
            for a in range(q-1):
                v[i][a] = math.log(fi[a])-(1/20)*log_sum

        w = np.zeros((len(x),len(x),q,q))
        obj = cls.from_parameters(v, w, **kwargs)
        obj.training_set = seq_file
        return obj


    @classmethod
    def from_parameters(cls, v, w, **kwargs): 
        """
            initialize MRF from pre-computed parameters
        """
        mrf = cls(v, w, **kwargs)
        if 'seq_file' in kwargs:
            mrf.training_set = kwargs['seq_file']
        if 'filename' in kwargs:
            if kwargs['filename'] is not None:
                mrf.to_msgpack(kwargs['filename'])
        return mrf


    def to_msgpack(self, filename=None):
        if filename is None:
            filename = self.name.replace('/','-')
        with open(str(filename), 'wb') as f:
            x_single = self.v[:,:20].reshape(self.ncol*20).tolist()
            x_pair = {}
            for i in range(self.ncol):
                for j in range(i + 1, self.ncol):
                    x_pair["{0}/{1}".format(i, j)] = {
                        "i": i,
                        "x": self.w[i, j, :, :].reshape(21 * 21).tolist(),
                        "j": j
                        }

                    out = {
                        "format": "ccm-1",
                        "ncol": self.ncol,
                        "x_single": x_single,
                        "x_pair": x_pair
                        }

            f.write(msgpack.packb(out, encoding="utf-8"))

        return pathlib.Path(filename)


    def get_w_norms(self):
        if hasattr(self, 'w_norms'):
            return self.w_norms
        else:
            w_norms = np.zeros((self.ncol, self.ncol))
            for i in range(0, self.ncol-1):
                for j in range(i+1, self.ncol):
                    w_norms[i][j] = self.get_w_norm_at_pos(i,j)
                    w_norms[j][i] = w_norms[i][j]
            self.w_norms = w_norms
            return w_norms



    def get_w_norm_at_pos(self, i, j):
        if hasattr(self, 'w_norms'):
            return self.w_norms[i][j]
        else:
            return math.sqrt(scalar_product(self.w[i][j],self.w[i][j]))


    def get_v_norms(self):
        return np.asarray([euclidean_norm(vi) for vi in self.v])


    def get_v_norm(self):
        return LA.norm(self.v)

    def get_w_norm(self):
        return LA.norm(self.w)


    def get_normalized_w_norms(self):
        w_norms = self.get_w_norms()
        m = np.mean(w_norms)
        s = np.std(w_norms)
        return (w_norms-m)/s

    def get_normalized_w_norm(self):
        return LA.norm(self.get_normalized_w_norms())


    def Hamiltonian(self, a):
        n = self.ncol
        s1=0
        for i in range(n-1):
            for j in range(i+1,n):
                s1+=self.w[i, j, code(a[i]), code(a[j])]

        s2=0
        for i in range(n):
            if a[i]!='-':
                s2+=self.v[i, code(a[i])]
        return s1+s2


    def Zi(self, x, i):
        """
            retourne la constante Zi à la position i pour la séquence x
        """
        n = self.ncol
        Zi = 0
        for a in range(21):
            sw=0
            for j in range(i+1,n):
                sw+=self.w[i, j, a, code(x[j])]
            Zi += math.exp(self.v[i, a] + sw)
        return Zi


    def log_pseudo_likelihood(self,x):
        """
            retourne le log du pseudo-likelihood de la séquence x pour le MRF
        """
        n = self.ncol
        p = 0
        for i in range(n-1):
            Zi = self.Zi(x, i)
            sw=0
            for j in range(i+1,n):
                sw+=self.w[i, j, code(x[i]), code(x[j])]
            p+=self.v[i,code(x[i])]+sw-math.log(Zi)
        return p
