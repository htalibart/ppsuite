""" ComPotts Object : contains MSA files locations, Potts models, etc. """

import os
import basic_modules.files_management as fm
from basic_modules.util import *
from basic_modules.tool_wrapper import *
from basic_modules.potts_model import *
from compotts.rescaling import *
from compotts.align_msas import *

class ComPotts_Object:

    @classmethod
    def from_hhblits_output(cls, seq_file, a3m_file, output_folder, input_folder=None, mrf_file=None, hhfilter_threshold=80, nb_sequences=1000, perform_trim=True, trimal_gt=0.8, trimal_cons=60, rescaling_function="identity", use_w=True, set_mrf=True, **kwargs):
        obj = cls()
        if 'name' in kwargs:
            obj.name = kwargs['name']
        else:
            obj.name = fm.get_name_from_first_sequence_name(seq_file)+"_hhblits"
        obj.seq_file = seq_file
        obj.a3m_file = a3m_file
        if input_folder is not None:
            obj.folder = input_folder
        else:
            if len(output_folder)>0:
                output_folder = fm.create_folder(output_folder)
            obj.folder = fm.create_folder(os.path.join(output_folder,obj.name))
        obj.aln_filtered = os.path.join(obj.folder,obj.name+"_filtered_80.a3m")
        if (not os.path.isfile(obj.aln_filtered)):
            call_hhfilter(obj.a3m_file, obj.aln_filtered, hhfilter_threshold)
        obj.reformat_file=os.path.join(obj.folder,obj.name+"_reformat.fasta")
        if (not os.path.isfile(obj.reformat_file)):
            call_reformat(obj.aln_filtered, obj.reformat_file)
        obj.aln_less = os.path.join(obj.folder,obj.name+"_less.fasta")
        if (not os.path.isfile(obj.aln_less)):
            fm.create_fasta_file_with_less_sequences(obj.reformat_file, obj.aln_less, nb_sequences)
        if perform_trim:
            print("trim ON")
            obj.colnumbering_file = os.path.join(obj.folder,obj.name+"_colnumbering.csv")
            obj.aln_trimmed = os.path.join(obj.folder,obj.name+"_trim_"+str(int(trimal_gt*100))+".fasta")
            if (not os.path.isfile(obj.aln_trimmed)):
                call_trimal(obj.aln_less, obj.aln_trimmed, trimal_gt, trimal_cons, obj.colnumbering_file)
            obj.trimal_ncol = fm.get_trimal_ncol(obj.colnumbering_file)
            obj.train_msa = obj.aln_trimmed
        else:
            print("trim OFF")
            obj.train_msa = obj.aln_less

        if set_mrf: # useful if we only want the MSA files
            if mrf_file is None:
                obj.mrf_file = os.path.join(obj.folder,obj.name+".mrf")
                if not os.path.isfile(obj.mrf_file):
                    obj.mrf = Potts_Model.from_training_set(obj.train_msa, obj.mrf_file, name=obj.name, **kwargs)
                else:
                    obj.mrf = Potts_Model.from_msgpack(obj.mrf_file, name=obj.name, **kwargs)
            else:
                obj.mrf_file = mrf_file
                obj.mrf = Potts_Model.from_msgpack(obj.mrf_file, name=obj.name, **kwargs)
            obj.original_mrf=obj.mrf
            if (rescaling_function!="identity"):
                print("rescaling MRF")
                obj.mrf = get_rescaled_mrf(obj.mrf, rescaling_function, use_w=use_w)
            else:
                print("using MRF as is (no rescaling)")
            obj.real_seq = fm.get_first_sequence_in_fasta_file(obj.seq_file).upper()
            obj.trimmed_seq = ''.join([obj.real_seq[col] for col in obj.get_real_positions([k for k in range(obj.mrf.ncol)])])
            return obj


    @classmethod
    def from_folder(cls, folder, **args):
         seq_file = fm.get_sequence_file_from_folder(folder)
         a3m_file = fm.get_a3m_file_from_folder(folder)
         mrf_file = fm.get_potts_model_file_from_folder(folder)
         return ComPotts_Object.from_hhblits_output(seq_file, a3m_file, folder, input_folder=folder, mrf_file=mrf_file, **args)



    @classmethod 
    def from_seq_file_to_one_hot(cls, seq_file, output_folder, rescaling_function="identity", use_w=True, **kwargs): # TODO tester
        obj = cls()
        if 'name' in kwargs:
            obj.name = kwargs['name']
        else:
            obj.name = fm.get_name_from_first_sequence_name(seq_file)+"_one_hot"
        obj.folder = fm.create_folder(os.path.join(output_folder,obj.name))
        obj.seq_file = seq_file
        obj.real_seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        obj.trimmed_seq = obj.real_seq
        obj.mrf_file = os.path.join(obj.folder,obj.name+".mrf")
        obj.mrf = Potts_Model.from_seq_file_to_one_hot(obj.seq_file, name=obj.name, **kwargs)
        obj.mrf.to_msgpack(obj.mrf_file)
        obj.original_mrf=obj.mrf
        if (rescaling_function!="identity"):
            print("rescaling MRF")
            obj.mrf = get_rescaled_mrf(obj.mrf, rescaling_function, use_w=use_w)
        else:
            print("using MRF as is (no rescaling)")
        obj.train_msa = obj.seq_file
        return obj

    @classmethod 
    def from_seq_file_with_submat(cls, seq_file, output_folder, rescaling_function="identity", use_w=True, **kwargs): # TODO tester
        obj = cls()
        if 'name' in kwargs:
            obj.name = kwargs['name']
        else:
            obj.name = fm.get_name_from_first_sequence_name(seq_file)+"_one_submat"
        obj.folder = fm.create_folder(os.path.join(output_folder,obj.name))
        obj.seq_file = seq_file
        obj.real_seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        obj.trimmed_seq = obj.real_seq
        obj.mrf_file = os.path.join(obj.folder,obj.name+".mrf")
        obj.mrf = Potts_Model.from_seq_file_with_submat(obj.seq_file, name=obj.name, **kwargs)
        obj.mrf.to_msgpack(obj.mrf_file)
        obj.original_mrf=obj.mrf
        if (rescaling_function!="identity"):
            print("rescaling MRF")
            obj.mrf = get_rescaled_mrf(obj.mrf, rescaling_function, use_w=use_w)
        else:
            print("using MRF as is (no rescaling)")
        obj.train_msa = obj.seq_file
        return obj


    @classmethod
    def from_seq_file_via_ccmpred(cls, seq_file, output_folder, rescaling_function = "identity", use_w=True, **kwargs):
        obj = cls()
        if 'name' in kwargs:    
            obj.name = kwargs['name']
        else:
            obj.name = fm.get_name_from_first_sequence_name(seq_file)+"_submat"
        obj.folder = fm.create_folder(os.path.join(output_folder,obj.name))
        obj.seq_file = seq_file
        obj.real_seq = fm.get_first_sequence_in_fasta_file(seq_file).upper()
        obj.trimmed_seq = obj.real_seq
        obj.mrf_file = os.path.join(obj.folder,obj.name+".mrf")
        obj.mrf = Potts_Model.from_training_set(obj.seq_file, obj.mrf_file, pc_submat="")
        obj.train_msa = obj.seq_file
        obj.original_mrf=obj.mrf
        if (rescaling_function!="identity"):
            print("rescaling MRF")
            obj.mrf = get_rescaled_mrf(obj.mrf, rescaling_function, use_w=use_w)
        else:
            print("using MRF as is (no rescaling)")

        return obj


    @classmethod
    def from_merge(cls, obj1, obj2, aligned_positions, output_folder, rescaling_function="identity", hhfilter_threshold=80, use_w=True, **kwargs):
        obj = cls()
        if "name" in kwargs:
            obj.name = kwargs['name']
        else:
            obj.name = '_'.join([obj1.name, obj2.name])
        obj.folder = fm.create_folder(os.path.join(output_folder,obj.name))
        obj.aln_unfiltered = os.path.join(obj.folder,obj.name+"_unfiltered.fasta")
        get_msas_aligned(aligned_positions, [obj1.train_msa, obj2.train_msa], obj.aln_unfiltered)
        obj.aln_filtered = os.path.join(obj.folder,obj.name+"_filtered.fasta")
        if (not os.path.isfile(obj.aln_filtered)):
            call_hhfilter(obj.aln_unfiltered, obj.aln_filtered, hhfilter_threshold)
        obj.train_msa = obj.aln_filtered
        obj.mrf_file = os.path.join(obj.folder,obj.name+".mrf")
        if os.path.isfile(obj.mrf_file):
            obj.mrf = Potts_Model.from_msgpack(obj.mrf_file, name=obj.name, **kwargs)
        else:
            obj.mrf = Potts_Model.from_training_set(obj.train_msa, obj.mrf_file, name=obj.name, **kwargs)
        obj.original_mrf=obj.mrf
        if (rescaling_function!="identity"):
            print("rescaling MRF")
            obj.mrf = get_rescaled_mrf(obj.mrf, rescaling_function, use_w=use_w)
        else:
            print("using MRF as is (no rescaling)")
        return obj


    def get_real_positions(self, positions_list): # TODO test
        if hasattr(self, 'trimal_ncol'):
            before_trimal = [self.trimal_ncol[j] for j in positions_list]
        else:
            before_trimal = positions_list
        if hasattr(self, 'reformat_file'):
            reformat_seq = fm.get_first_sequence_in_fasta_file(self.reformat_file)
            small_to_real_list = get_small_to_real_list(self.real_seq, reformat_seq)
        else:
            small_to_real_list = [k for k in range(len(self.real_seq))]
        real_positions = [small_to_real_list[j] for j in before_trimal]
        return real_positions


    def get_real_letter_at_trimmed_pos(self, position): # TODO check if deprecated
        return self.trimmed_seq[position]

