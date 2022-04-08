import os
import pathlib
import numpy as np
import numpy.linalg as LA
import math
import msgpack
import json
import subprocess

import ccmpred.weighting
import ccmpred.io
import ccmpred.counts

from comutils.util import *
from comutils import files_management as fm
from comutils import pseudocounts
from comutils.adabmdca_to_ccmpredpy import *

from rmfdca import __main__ as rmfdca_main

POSSIBLE_CCMPRED_OPTIONS = ["wt-simple", "wt-simple", "wt-uniform", "wt-cutoff", "reg-lambda-single", "reg-lambda-pair-factor", "reg-L2", "reg-noscaling", "reg-scale-by-L", "v-center", "v-zero", "max-gap-pos", "max-gap_seq", "pc-uniform", "pc-submat", "pc-constant", "pc-none", "pc-single-count", "pc-pair-count", "maxit", "ofn-pll", "ofn-cd", "pc-pair-submat", "persistent", "no-decay", "nr-markov-chains"]

MFDCA_ALPHABET="ACDEFGHIKLMNPQRSTVWY-"


class Potts_Model:

    def __init__(self, v, w, **kwargs):
        self.v = v.astype(np.float32)
        self.ncol = len(self.v)
        self.w = w.astype(np.float32)
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = "Billy"


    @classmethod
    def from_msgpack(cls, binary_file, **kwargs):
        """
            initialize MRF from msgpack file
        """
        with open(str(binary_file), 'rb') as data_file:
            df = msgpack.unpackb(data_file.read())
            """
            CCMpredPy: df is a dictionary
                b'format'
                b'ncol'
                b'x_single': np.array(ncol, 20)) 
                b'x_pair' : np.array(ncol, ncol, 21, 21)
                b'meta'
            """
        print("getting Potts model from "+str(binary_file))
        fm.check_if_file_ok(binary_file)
        ncol = df[b'ncol']
        q = int(len(df[b'x_single'])/ncol) # real q is considered to be v.shape[1]
        v = np.array(df[b'x_single']).reshape((ncol,q))
        w = np.zeros((ncol, ncol, q, q))
        for p in df[b'x_pair'].values():
            i = p[b'i']
            j = p[b'j']
            q_w = int(np.sqrt(len(p[b'x'])))# CCMpredPy gives 21x21 w matrices even though q=20...
            mat = np.array(p[b'x']).reshape((q_w,q_w))
            w[i, j, :, :] = mat[:q,:q]
            w[j, i, :, :] = mat.T[:q,:q]
        if 'name' not in kwargs:
            kwargs['name'] = str(binary_file).replace('/','-')
            if kwargs['name'].startswith('-'):
                kwargs['name'] = kwargs['name'][1:]
        mrf = cls(v, w, **kwargs)
        mrf.binary_file = binary_file
        return mrf



    @classmethod
    def from_adabmdca_file(cls, adabmdca_file, binary_file=None, **kwargs):
        with open(str(adabmdca_file), 'r') as adaf:
            csv_reader = csv.reader(adaf, delimiter=' ')
            J_lines = []
            h_lines = []
            for line in csv_reader:
                if line[0]=='J':
                    J_lines.append([int(ind) for ind in line[1:5]]+[float(line[5])])
                else:
                    h_lines.append([int(ind) for ind in line[1:3]]+[float(line[3])])
            ncol = len(h_lines)//21
            v = np.zeros((ncol, 21))
            for h_line in h_lines:
                v[h_line[0]][h_line[1]] = h_line[2]
            w = np.zeros((ncol, ncol, 21, 21))
            for J_line in J_lines:
                w[J_line[0],J_line[1],J_line[2],J_line[3]]=J_line[4]
        mrf = cls.from_parameters(v, w, **kwargs)
        if binary_file is not None:
            mrf.to_msgpack(binary_file)
        return mrf



    @classmethod
    def from_training_set(cls, aln_file, binary_file, write_readme=True, readme_file=None, inference_method='CCMpredPy',  wt_cutoff=0.8, lattice=False, reg_lambda_w_mfdca=1, pc_tau_mfdca=0.5, shrinkage_coeff=0.5, mfdca_pc=False, **kwargs):
        """
            initialize Potts model from train MSA file
        """
        fm.check_if_file_ok(aln_file)

        if inference_method=='CCMpredPy':
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

        elif inference_method=='adabmDCA':
            adabmDCA_output_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
            adabmDCA_output_folder.mkdir()
            #call = "adabmDCA -f "+str(aln_file)+" -z 1 -m 500 -c 1e-2"
            call = "adabmDCA -f "+str(aln_file)+" -z 1 -m 500 -c 1e-2 -i 1"
            subprocess.Popen(call, shell=True, cwd=str(adabmDCA_output_folder)).wait()
            weights_file = adabmDCA_output_folder/'Parameters_conv_zerosum_nolabel.dat'
            if not os.path.exists(weights_file):
                raise Exception("No adabmDCA output file "+str(weights_file))
            mrf = cls.from_adabmdca_file(weights_file, binary_file=binary_file)
            shutil.rmtree(adabmDCA_output_folder)


        elif inference_method=='mfDCA':
            v, w = rmfdca_main.infer_parameters_for_msa(aln_file, reg_lambda_w=reg_lambda_w_mfdca, pc_tau=pc_tau_mfdca, shrinkage_coeff=shrinkage_coeff, lattice=lattice, mfdca_pc=mfdca_pc) 
            mrf = cls.from_parameters(v, w)
            mrf.to_msgpack(binary_file)
            

        else:
            raise Exception("Cannot infer Potts models with "+inference_method+" (available: CCMpredPy, adabmDCA, mfDCA)")

        mrf.training_set = pathlib.Path(aln_file)
        return mrf


    @classmethod
    def from_sequence_to_one_hot(cls, seq, seq_file=None, q=21, **kwargs):
        """ one hot encoding """
        x = code_whole_seq(seq)
        v = np.zeros((len(x),q))
        for i in range(len(x)):
            v[i,x[i]]=1
        w = np.zeros((len(x),len(x),q,q))
        for i in range(len(x)):
            for j in range(len(x)):
                w[i,j,x[i],x[j]] = 1
        obj = cls.from_parameters(v, w, **kwargs)
        obj.training_set = seq_file
        return obj

    @classmethod
    def from_sequence_file_to_one_hot(cls, seq_file, q=21, **kwargs):
        """ one hot encoding """
        fm.check_if_file_ok(seq_file)
        seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        return cls.from_sequence_to_one_hot(seq, seq_file=seq_file, q=q, **kwargs)


    @classmethod
    def from_sequence_with_submat_freq(cls, seq, seq_file=None, tau=0.5, q=21, **kwargs):
        """ substitution matrix pseudocounts to frequencies"""
        x = code_whole_seq(seq)
        
        v = np.zeros((len(x),q))
        for i in range(len(x)):
            for a in range(q-1):
                v[i,a] = (1-tau)*(a==x[i]) + tau*pseudocounts.get_cond_proba(a,x[i])

        w = np.zeros((len(x),len(x),q,q))
        for i in range(len(x)):
            for j in range(len(x)):
                w[i,j,x[i],x[j]] = 1
        obj = cls.from_parameters(v, w, **kwargs)
        obj.training_set = seq_file
        return obj

    @classmethod
    def from_sequence_file_with_submat_freq(cls, seq_file, tau=0.5, q=21, **kwargs):
        """ substitution matrix pseudocounts to frequencies"""
        fm.check_if_file_ok(seq_file)
        seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        return cls.from_sequence_with_submat_freq(seq, seq_file, tau=tau, q=q, **kwargs)


    @classmethod
    def from_sequence_with_submat(cls, seq, seq_file=None, tau=0.5, q=21, **kwargs): #TODO check
        """ substitution matrix pseudocounts """
        x = code_whole_seq(seq)
        v = np.zeros((len(x),q))

        for i in range(len(x)):
            fi = np.zeros(q-1)
            log_sum=0
            for a in range(q-1):
                fi[a] = (1-tau)*(a==x[i]) + tau*pseudocounts.get_cond_proba(a,x[i])
                log_sum+=math.log(fi[a])
            for a in range(q-1):
                v[i][a] = math.log(fi[a])-(1/20)*log_sum

        w = np.zeros((len(x),len(x),q,q))
        for i in range(len(x)):
            for j in range(len(x)):
                w[i,j,x[i],x[j]] = 1
        obj = cls.from_parameters(v, w, **kwargs)
        obj.training_set = seq_file
        return obj


    
    @classmethod
    def from_sequence_file_with_submat(cls, seq_file, tau=0.5, **kwargs):
        """ substitution matrix pseudocounts """
        fm.check_if_file_ok(seq_file)
        seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        return cls.from_sequence_with_submat(seq, seq_file, **kwargs)


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


    @classmethod
    def zero_fill(cls, length, q=21, **kwargs):
        return cls.from_parameters(v=np.zeros((length,q)), w=np.zeros((length,length,q,q)), **kwargs)


    def to_msgpack(self, filename=None):
        if filename is None:
            filename = self.name.replace('/','-')
        q = self.v.shape[1]
        with open(str(filename), 'wb') as f:
            x_single = self.v.reshape(self.ncol*q).tolist()
            x_pair = {}
            for i in range(self.ncol):
                for j in range(i + 1, self.ncol):
                    x_pair["{0}/{1}".format(i, j)] = {
                        "i": i,
                        "x": self.w[i, j, :, :].reshape(q * q).tolist(),
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
        w_norms = np.zeros((self.ncol, self.ncol))
        for i in range(0, self.ncol):
            for j in range(0, self.ncol):
                w_norms[i][j] = self.get_w_norm_at_pos(i,j)
                #w_norms[j][i] = w_norms[i][j]
        return w_norms



    def get_w_norm_at_pos(self, i, j):
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
            if (a[i]!='-') or (self.w.shape[2]==21):
                for j in range(i+1,n):
                    if (a[j]!='-') or (self.w.shape[2]==21):
                        s1+=self.w[i, j, code(a[i]), code(a[j])]

        s2=0
        for i in range(n):
            if (a[i]!='-') or (self.v.shape[1]==21):
                s2+=self.v[i, code(a[i])]
        return s1+s2


    def Zi(self, x, i):
        """ pseudo-likelihood constant Zi at position i for sequence x """
        n = self.ncol
        Zi = 0
        q = self.v.shape[1]
        for a in range(q):
            sw=0
            for j in range(i+1,n):
                sw+=self.w[i, j, a, code(x[j])]
            Zi += math.exp(self.v[i, a] + sw)
        return Zi


    def log_pseudo_likelihood(self,x):
        """ log pseudo-likelihood of sequence @x """
        n = self.ncol
        p = 0
        for i in range(n):
            Zi = self.Zi(x, i)
            sw=0
            for j in range(i+1,n):
                if j<n:
                    sw+=self.w[i, j, code(x[i]), code(x[j])]
            p+=self.v[i,code(x[i])]+sw-math.log(Zi)
        return p


    def insert_null_position_at(self, pos, v_null=None):
        """ insert null column at position @pos where the field vector is @v_null and all w coupled with @pos are 0 """
        q = self.v.shape[1]
        if v_null is None:
            v_null = np.zeros((1,q))
        self.ncol = self.ncol+1
        self.v = np.concatenate((self.v[:pos],v_null,self.v[pos:]))
        new_w = np.zeros((self.ncol,self.ncol,q,q))
        for i in range(self.ncol-1):
            for j in range(i+1,self.ncol):
                if (i==pos) or (j==pos):
                    new_w[i,j]=np.zeros((q,q))
                else:
                    new_w[i,j]=self.w[i-(i>pos),j-(j>pos)]
                new_w[j,i] = new_w[i,j]
        self.w = new_w


    def insert_null_positions_to_complete_mrf_pos(self, mrf_pos_to_seq_pos, sequence_length, v_null=None):
        """ insert null column at each position in the sequence which is not in the Potts model """
        q = self.v.shape[1]
        if v_null is None:
            v_null = np.zeros((1,q))
        for pos_in_seq in range(sequence_length):
            if not pos_in_seq in mrf_pos_to_seq_pos:
                self.insert_null_position_at(pos_in_seq, v_null)

    def insert_vi_to_complete_mrf_pos(self, mrf_pos_to_seq_pos, sequence_length, v_fill):
        """ insert appropriate column of v_fill at each position in the sequence which is not in the Potts model """
        for pos_in_seq in range(sequence_length):
            if not pos_in_seq in mrf_pos_to_seq_pos:
                v_i = v_fill[pos_in_seq].reshape((1,q))
                self.insert_null_position_at(pos_in_seq, v_null=v_i)


    def insert_vi_star_gapped_to_complete_mrf_pos(self, mrf_pos_to_seq_pos, sequence_length, msa_file_before_trim, nb_pc_for_v_star=1, wt_cutoff=0.8):
        """ insert v* column at each position in the sequence which is not in the Potts model """
        v_star = compute_v_star(msa_file_before_trim, wt_cutoff, nb_pc_for_v_star)
        q = self.v.shape[1]
        if (q==21):
            v_star[:, 20]=0
        self.insert_vi_to_complete_mrf_pos(mrf_pos_to_seq_pos, sequence_length, v_star)


    def insert_vi_with_blosum_pseudocounts_to_complete_mrf_pos(self, mrf_pos_to_seq_pos, sequence_length, msa_file_before_trim, freq_gap_min, pc_tau):
        v_fill = compute_v_with_blosum_pseudocounts_for_gaps(msa_file_before_trim, freq_gap_min, pc_tau)
        self.insert_vi_to_complete_mrf_pos(mrf_pos_to_seq_pos, sequence_length, v_fill)



    def change_gauge_l2_zero_to_l2_center(self, v_star, reg_lambda_single=None, reg_lambda_pair_factor=None):
        q = self.v.shape[1]
        if (reg_lambda_single is None) or (reg_lambda_pair_factor is None):
            reg_ratio = np.sum(self.w[0,1,0])/self.v[0,0]
        else:
            reg_ratio = reg_lambda_single/(reg_lambda_pair_factor*(self.ncol-1))
        new_w = np.zeros_like(self.w)
        for i in range(self.ncol-1):
            for j in range(i+1,self.ncol):
                for a in range(q):
                    for b in range(q):
                        new_w[i,j,a,b] = self.w[i,j,a,b] - (reg_ratio/q)*(v_star[i,a] + v_star[j,b])
                        new_w[j,i,b,a] = new_w[i,j,a,b]
        new_v = self.v + ((self.ncol-1)*reg_ratio/q)*v_star
        return Potts_Model.from_parameters(new_v, new_w)


    def change_gauge_l2_center_to_l2_zero(self, v_star, reg_lambda_single=None, reg_lambda_pair_factor=None):
        q = self.v.shape[1]
        if (reg_lambda_single is None) or (reg_lambda_pair_factor is None):
            reg_ratio = np.sum(self.w[0,1,0])/(self.v[0,0]-v_star[0,0])
        else:
            reg_ratio = reg_lambda_single/(reg_lambda_pair_factor*(self.ncol-1))
        new_w = np.zeros_like(self.w)
        for i in range(self.ncol-1):
            for j in range(i+1,self.ncol):
                for a in range(q):
                    for b in range(q):
                        new_w[i,j,a,b] = self.w[i,j,a,b] + (reg_ratio/q)*(v_star[i,a] + v_star[j,b])
                        new_w[j,i,b,a] = new_w[i,j,a,b]
        new_v = self.v - ((self.ncol-1)*reg_ratio/q)*v_star
        return Potts_Model.from_parameters(new_v, new_w)


    def change_gauge_zero_sum_to_l2_center(self, v_star, lv=1, lw=1):
        q = self.v.shape[1]
        reg_ratio = lv/lw
        new_w = np.zeros_like(self.w)
        for i in range(self.ncol-1):
            for j in range(i+1,self.ncol):
                for a in range(q):
                    for b in range(q):
                        new_w[i,j,a,b] = self.w[i,j,a,b] + (reg_ratio/q)*(self.v[i,a]-v_star[i,a] + self.v[j,b]-v_star[j,b])
                        new_w[j,i,b,a] = new_w[i,j,a,b]
        new_v = self.v - ((self.ncol-1)*reg_ratio/q)*(self.v-v_star)
        return Potts_Model.from_parameters(new_v, new_w)


    def apply_zero_sum_gauge(self):
        v = self.v
        w = self.w
        zv = np.zeros_like(v)
        zw = np.zeros_like(w)
        L = v.shape[0]
        q = v.shape[1]
        average_v = np.mean(v, axis=1)
        average_w = np.mean(w, axis=(2,3))
        average_w_a = np.mean(w, axis=2)
        average_w_b = np.mean(w, axis=3)
        for i in range(L):
            for a in range(q):
                zv[i,a] = v[i,a]-average_v[i]+np.sum([average_w_b[i,j,a]-average_w[i,j] for j in range(L) if j!=i])
        for i in range(L):
            for j in range(L):
                for a in range(q):
                    for b in range(q):
                        zw[i,j,a,b] = w[i,j,a,b]-average_w_b[i,j,a]-average_w_a[i,j,b]+average_w[i,j]
                #zw[j,i] = np.transpose(zw[i,j])
        return Potts_Model.from_parameters(zv, zw)
