""" ComPotts Object : contains MSA files locations, Potts models, etc. """

import os
import time
import basic_modules.files_management as fm
from basic_modules.util import *
from basic_modules.tool_wrapper import *
from basic_modules.potts_model import *
from compotts.rescaling import *
from compotts.align_msas import *

class ComPotts_Object:

    def __init__(self, mrf=None, mrf_file=None, name=None, seq_file=None, aln_fasta=None, a3m_file=None, input_folder=None, nb_sequences=1000, use_less_sequences=True, hhfilter_threshold=80, perform_filter=True, trimal_gt=0.8, trimal_cons=60, pc_count_factor=1000, reg_lambda_pair_factor=30, trim_alignment=True, rescaling_function="identity", use_w=True, mrf_type=None):

        # SEQ_FILE
        if seq_file is not None:
            self.seq_file = seq_file
        elif input_folder is not None:
            self.seq_file = fm.get_sequence_file_from_folder(input_folder)
        else:
            self.seq_file = None

        # EXISTING ALIGNMENT FASTA FORMAT
        self.aln_fasta = aln_fasta

        # EXISTING A3M FILE
        if a3m_file is not None:
            self.a3m_file = a3m_file
        elif input_folder is not None:
            self.a3m_file = fm.get_a3m_file_from_folder(input_folder)
        else:
            self.a3m_file = None

        # NAME
        if name is not None:
            self.name = name
        elif self.seq_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.seq_file)
        elif self.aln_fasta is not None:
            self.name = fm.get_first_sequence_clean_name(self.aln_fasta)
        elif self.a3m_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.a3m_file)
        elif mrf is not None:
            self.name = mrf.name
        else:
            self.name = "Billy_"+time.strftime("%Y%m%d-%H%M%S")


        # FOLDER
        if input_folder is not None:
            self.folder = input_folder
        else:
            self.folder = fm.create_folder(self.name)

 
        # REFORMAT A3M_FILE
        if (self.aln_fasta is None) and (self.a3m_file is not None):
            self.aln_file = self.folder/(name+"_reformat.fasta")
            call_reformat(self.a3m_file, self.aln_file)

        # SEQUENCE
        if self.seq_file is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.seq_file)
        elif self.aln_fasta is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.aln_fasta)


        # FILTER
        if (self.aln_fasta is not None) and (perform_filter):
            old_aln_file = self.aln_file
            self.aln_file = self.folder/(name+"_filter_"+str(hhfilter_threshold)+".fasta")
            if not self.aln_file.is_file():
                call_hhfilter(old_aln_file, self.aln_file, hhfilter_threshold):


        # USE LESS SEQUENCES
        if (self.aln_fasta is not None) and (use_less_sequences):
            old_aln_file = self.aln_file
            self.aln_file = self.folder/(name+"_less.fasta")
            if not self.aln_file.is_file():
                fm.create_fasta_file_with_less_sequences(old_aln_file, self.aln_file, nb_sequences)


        # TRIM ALIGNMENT
        if (self.aln_fasta is not None) and (trim_alignment):
            old_aln_file = self.aln_file
            self.aln_file = self.folder/(name+"_trim_"+str(trimal_gt*100)+".fasta"
            colnumbering_file = self.folder/(name+"_colnumbering.csv")
            if not self.aln_file.is_file():
                call_trimal(old_aln_file, self.aln_file, trimal_gt, cons, colnumbering_file)
            self.real_aln_pos = fm.get_trimal_ncol(colnumbering_file)
        elif (self.aln_fasta is not None):
            nb_pos = fm.get_nb_columns_in_alignment(self.aln_fasta) 
            self.real_aln_pos = [pos for pos in range(nb_pos)]


        # ALIGNMENT POSITIONS -> SEQUENCE POSITIONS
        if (self.aln_fasta is not None) and (self.sequence is not None):
            aln_first_seq = fm.get_first_sequence_in_fasta_file(self.aln_fasta)
            seq_aln_pos = fm.get_small_to_real_list(self.sequence, aln_first_seq)
            self.real_seq_pos = seq_to_aln_pos[pos for pos in self.real_aln_pos]
        elif (self.sequence is not None):
            self.real_seq_pos = [pos for pos in range(len(self.sequence))]


        # DIFFÉRENTES FAÇONS DE CRÉER UN MRF
        if mrf_type is not None:
            self.mrf_type=mrf_type
        elif self.aln_fasta is not None: # if aln_fasta exists, MRF is trained in a standard way
            self.mrf_type="standard"
        elif self.seq_file is not None: # if aln_fasta doesn't exist but we have a sequence file, it is used to train the MRF
           self.mrf_type="one_submat"


        # TRAINING SET EN FONCTION DES DIFFÉRENTES FAÇONS DE CRÉER UN MRF
        if (self.mrf_type=="standard"):
            if (self.aln_fasta) is not None:
                self.training_set = self.aln_fasta
            else:
                print("Missing alignment file !")
        elif (self.mrf_type=="one_submat") or (self.mrf_type=="one_hot"):
            if self.seq_file is not None:
                self.training_set = self.seq_file
            elif self.sequence is not None:
                self.seq_file = output_folder/(self.name+".fasta")
                create_seq_fasta(self.sequence, self.seq_file, seq_name=self.name)
                self.training_set = self.seq_file
            else:
                print("Missing sequence file !")

        # MRF
        if mrf is not None:
            self.mrf = mrf
        else:
            if mrf_file is not None:
                self.mrf_file = mrf_file
                self.mrf = Potts_Model.from_msgpack(mrf_file, **kwargs)
            else:
                self.mrf_file = (self.folder)/(name+"_"+mrf_type".fasta")
                if (self.mrf_type=="standard"):
                    self.mrf = Potts_Model.from_training_set(self.training_set, self.mrf_file, **kwargs)
                elif (self.mrf_type=="one_hot"):
                    self.mrf = Potts_Model.from_seq_file_to_one_hot(self.training_set, **kwargs)
                elif (self.mrf_type=="one_submat"):
                    self.mrf = Potts_Model.from_seq_file_with_submat(self.training_set, **kwargs)
        self.original_mrf = mrf
        if (rescaling_function!="identity"):
            print("rescaling MRF")
            self.mrf = get_rescaled_mrf(self.mrf, rescaling_function, use_w=use_w)



    @classmethod
    def from_merge(cls, obj1, obj2, aligned_positions, output_folder, **kwargs):
        if "name" in kwargs:
            name = kwargs['name']
        else:
            name = '_'.join([obj1.name, obj2.name])

        if input_folder is None:
            input_folder = fm.create_folder(name)

        aln_unfiltered = input_folder/(name+"_unfiltered.fasta")
        if not aln_unfiltered.is_file():
            get_msas_aligned(aligned_positions, [obj1.training_set, obj2.training_set], aln_unfiltered)
        obj = cls(self, name=name, aln_fasta=aln_unfiltered, input_folder=input_folder, **kwargs)

        return obj
