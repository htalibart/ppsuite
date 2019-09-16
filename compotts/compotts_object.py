""" ComPotts Object : contains MSA files locations, Potts models, etc. """

import os
import time
import basic_modules.files_management as fm
from basic_modules.util import *
from basic_modules.tool_wrapper import *
from basic_modules.potts_model import *
from compotts.rescaling import *
from compotts.align_msas import *
from compotts.manage_positions import *

class ComPotts_Object:

    def __init__(self, mrf=None, potts_model_file=None, name=None, sequence_file=None, aln_fasta=None, a3m_file=None, input_folder=None, nb_sequences=1000, use_less_sequences=True, hhfilter_threshold=80, perform_filter=True, trimal_gt=0.8, trimal_cons=60, pc_count=1000, reg_lambda_pair_factor=30, trim_alignment=True, rescaling_function="identity", use_w=True, mrf_type=None, hhblits_database=None, **kwargs):

        self.folder = input_folder

        # MRF
        self.potts_model_file = potts_model_file
        if (self.potts_model_file is None) and (input_folder is not None):
            self.potts_model_file = fm.get_potts_model_file_from_folder(input_folder)

        # SEQ_FILE
        if sequence_file is not None:
            self.sequence_file = sequence_file
        elif input_folder is not None:
            self.sequence_file = fm.get_sequence_file_from_folder(input_folder)
        else:
            self.sequence_file = None

        # EXISTING ALIGNMENT FASTA FORMAT
        self.aln_fasta = aln_fasta
        if self.aln_fasta is None:
            self.aln_fasta = fm.get_existing_training_set(input_folder, trimal_gt)

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
        elif self.sequence_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.sequence_file)
        elif self.aln_fasta is not None:
            self.name = fm.get_first_sequence_clean_name(self.aln_fasta)
        elif self.a3m_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.a3m_file)
        elif mrf is not None:
            self.name = mrf.name
        elif input_folder is not None:
            self.name = input_folder.stem
        else:
            self.name = "Billy_"+time.strftime("%Y%m%d-%H%M%S")



        # IF WE DON'T HAVE AN A3M FILE AND WE WANT ONE
        if (self.a3m_file is None) and (self.aln_fasta is None) and (mrf_type=="standard") and (self.potts_model_file is None):
            self.a3m_file = self.get_folder()/(self.name+".a3m")
            if hhblits_database is None:
                raise Exception("Must specify a database for hhblits (option -d)")
            else:
                call_hhblits(self.sequence_file, self.a3m_file, hhblits_database, **kwargs)

 
        # REFORMAT A3M_FILE
        if (self.aln_fasta is None) and (self.a3m_file is not None):
            self.a3m_reformat = self.get_folder()/(self.name+"_reformat.fasta")
            if (not self.a3m_reformat.is_file()) and (self.potts_model_file is None):
                call_reformat(self.a3m_file, self.a3m_reformat)
            self.aln_fasta = self.a3m_reformat


        # FILTER
        if (self.aln_fasta is not None) and (perform_filter) and (self.potts_model_file is None):
            old_aln_fasta = self.aln_fasta
            self.aln_fasta = self.get_folder()/(self.name+"_filtered_"+str(hhfilter_threshold)+".fasta")
            if (not self.aln_fasta.is_file()):
                call_hhfilter(old_aln_fasta, self.aln_fasta, hhfilter_threshold)


        # USE LESS SEQUENCES
        if (self.aln_fasta is not None) and (use_less_sequences) and (self.potts_model_file is None):
            old_aln_fasta = self.aln_fasta
            self.aln_fasta = self.get_folder()/(self.name+"_less.fasta")
            if (not self.aln_fasta.is_file()):
                fm.create_fasta_file_with_less_sequences(old_aln_fasta, self.aln_fasta, nb_sequences)


        # TRIM ALIGNMENT
        if (self.aln_fasta is not None) and (trim_alignment):
            colnumbering_file = self.get_folder()/(self.name+"_colnumbering.csv")
            if (self.potts_model_file is None):
                old_aln_fasta = self.aln_fasta
                self.aln_fasta = self.get_folder()/(self.name+"_trim_"+str(int(trimal_gt*100))+".fasta")
                if (not self.aln_fasta.is_file()):
                    call_trimal(old_aln_fasta, self.aln_fasta, trimal_gt, trimal_cons, colnumbering_file)
            self.real_aln_pos = fm.get_trimal_ncol(colnumbering_file)
        elif (self.aln_fasta is not None):
            nb_pos = fm.get_nb_columns_in_alignment(self.aln_fasta) 
            self.real_aln_pos = [pos for pos in range(nb_pos)]


        # SEQUENCE
        if self.sequence_file is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.sequence_file)
        elif self.aln_fasta is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.aln_fasta)
        else:
            self.sequence = None



        # ALIGNMENT POSITIONS -> SEQUENCE POSITIONS
        if (self.aln_fasta is not None) and (self.sequence is not None):
            if (hasattr(self, "a3m_reformat")):
                aln_first_seq = fm.get_first_sequence_in_fasta_file(self.a3m_reformat)
                seq_aln_pos = get_real_pos_list(self.sequence, aln_first_seq)
                self.real_seq_pos = [seq_aln_pos[pos] for pos in self.real_aln_pos]
            else:
                self.real_seq_pos = self.real_aln_pos
        elif (self.sequence is not None):
            self.real_seq_pos = [pos for pos in range(len(self.sequence))]
            self.real_aln_pos = self.real_seq_pos


        # DIFFÉRENTES FAÇONS DE CRÉER UN MRF
        if mrf_type is not None:
            self.mrf_type=mrf_type
        elif self.aln_fasta is not None: # if aln_fasta exists, MRF is trained in a standard way
            self.mrf_type="standard"
        elif self.sequence_file is not None: # if aln_fasta doesn't exist but we have a sequence file, it is used to train the MRF
           self.mrf_type="one_submat"
        else:
            self.mrf_type=mrf_type


        # TRAINING SET EN FONCTION DES DIFFÉRENTES FAÇONS DE CRÉER UN MRF
        self.training_set = None
        if (self.mrf_type=="standard"):
            if (self.aln_fasta) is not None:
                self.training_set = self.aln_fasta
            else:
                print("Missing alignment file !")
        elif (self.mrf_type=="one_submat") or (self.mrf_type=="one_hot"):
            if self.sequence_file is not None:
                self.training_set = self.sequence_file
            elif self.sequence is not None:
                self.sequence_file = self.folder/(self.name+".fasta")
                fm.create_seq_fasta(self.sequence, self.sequence_file, seq_name=self.name)
                self.training_set = self.sequence_file
            else:
                print("Missing sequence file !")


        # MRF 
        if mrf is not None:
            self.mrf = mrf
        else:
            if self.potts_model_file is not None:
                self.mrf = Potts_Model.from_msgpack(self.potts_model_file, **kwargs)
            else:
                if (self.training_set is not None):
                    self.potts_model_file = (self.get_folder())/(self.name+"_"+self.mrf_type+".mrf")
                    if (self.mrf_type=="standard"):
                        self.mrf = Potts_Model.from_training_set(self.training_set, self.potts_model_file, pc_count=pc_count, reg_lambda_pair_factor=reg_lambda_pair_factor, **kwargs)
                    elif (self.mrf_type=="one_hot"):
                        self.mrf = Potts_Model.from_sequence_file_to_one_hot(self.training_set, **kwargs)
                    elif (self.mrf_type=="one_submat"):
                        self.mrf = Potts_Model.from_sequence_file_with_submat(self.training_set, **kwargs)
                    os.system("cp "+str(self.training_set)+" "+str(self.get_folder()/(self.name+"_training_set.fasta")))
                else:
                    raise Exception("Need a training set")
            self.original_mrf = self.mrf
        if (rescaling_function!="identity"):
            print("rescaling MRF")
            self.mrf = get_rescaled_mrf(self.mrf, rescaling_function, use_w=use_w)



    @classmethod
    def from_merge(cls, obj1, obj2, aligned_positions, output_folder, **kwargs):
        """ Creating a ComPotts object by merging two ComPotts objects """
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


    def get_folder(self):
        if not hasattr(self, "folder"):
            self.folder = fm.create_folder(self.name)
        elif self.folder is None:
            self.folder = fm.create_folder(self.name)
        if (not self.folder.is_dir()):
            self.folder = fm.create_folder(self.folder)
        return self.folder
        

    def get_seq_positions(self, positions):
        return [self.real_seq_pos[pos] for pos in positions]
