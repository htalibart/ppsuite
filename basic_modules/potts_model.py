import os
import pathlib
import numpy as np
import math
import msgpack

from basic_modules.util import *
from basic_modules import files_management as fm
from basic_modules import pseudocounts

POSSIBLE_CCMPRED_OPTIONS = ["wt-simple", "wt-simple", "wt-uniform", "wt-cutoff", "reg-lambda-single", "reg-lambda-pair-factor", "reg-L2", "reg-noscaling", "reg-scale-by-L", "v-center", "v-zero", "max-gap-pos", "max-gap_seq", "pc-uniform", "pc-submat", "pc-constant", "pc-none", "pc-count", "pc-pair-count", "maxit"] # TODO mettre toutes les options

class Potts_Model:

    def __init__(self, v, w, **kwargs):
        self.v = v
        self.ncol = len(self.v)
        self.w = w
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = "Billy"


    @classmethod
    def from_msgpack(cls, binary_file, **kwargs):
        """
            initialize MRF from msgpack file given by CCMpredPy
        """
        with open(binary_file, 'rb') as data_file:
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
        mrf = cls(v, w, **kwargs)
        mrf.binary_file = binary_file
        return mrf

    @classmethod
    def from_training_set(cls, aln_file, binary_file, **kwargs):
        """
            initialize MRF from training set
        """
        call = "ccmpred "+str(aln_file)+ " -b "+str(binary_file)
        for key_arg in kwargs:
            arg_ccm = key_arg.replace('_', '-')
            if arg_ccm in POSSIBLE_CCMPRED_OPTIONS:
                call+=" --"+arg_ccm+" "+str(kwargs[key_arg])
        os.system(call)
        if not os.path.exists(binary_file):
            raise Exception("CCMpredPy wasn't able to infer the MRF. Protein is probably too long ?")
        mrf = cls.from_msgpack(binary_file)
        mrf.training_set = pathlib.Path(aln_file)
        return mrf


    @classmethod
    def from_seq_file_to_one_hot(cls, seq_file, **kwargs):
        """ one hot encoding """
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
    def from_seq_file_with_submat(cls, seq_file, npc=1, **kwargs): # TODO checker pourquoi le MRF est différent en inférant avec CCMpredPy
        """ substitution matrix pseudocounts """
        seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        tau = npc/(1+npc)
        x = code_whole_seq(seq)
        v = np.zeros((len(x), q))
        v = np.zeros((len(x),q))
        p_submat = pseudocounts.get_probas()
        for i in range(len(x)):
            sum_b = (1/(q-1))*sum([math.log((1-tau)*(x[i]==b)+tau*p_submat[b]) for b in range(q-1)])
            for a in range(q-1):
                v[i,a] = math.log((1-tau)*(x[i]==a)+tau*p_submat[a])-sum_b
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
        return mrf


    def to_msgpack(self, filename=None):
        if filename is None:
            filename = self.name.replace('/','-')
        with open(filename, 'wb') as f:
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
        return [euclidean_norm(vi) for vi in self.v]


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


    def get_rescaled_version(self, v_rescaling_function, w_rescaling_function):
        q=len(ALPHABET)
        t_v = np.zeros((mrf.ncol, q))
        for i in range(mrf.ncol):
            for a in range(q):
                t_v[i][a] = v_rescaling_function(self.v[i][a])
        t_w = np.zeros((mrf.ncol, mrf.ncol, q, q))
        for i in range(mrf.ncol):
                for j in range(mrf.ncol):
                    for a in range(q):
                        for b in range(q):
                            t_w_[i][j][a][b] = w_rescaling_function(self.w[i][j][a][b])
        t_mrf = Potts_Model.from_parameters(self.ncol, t_v, t_w, name=self.name+"_rescaled")
        return t_mrf
