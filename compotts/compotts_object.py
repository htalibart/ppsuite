""" ComPotts Object : contains MSA files locations, Potts models, etc. """

import os
import time
import argparse
import sys
import basic_modules.files_management as fm
from basic_modules.util import *
from basic_modules.tool_wrapper import *
from basic_modules.potts_model import *
from compotts.rescaling import *
from compotts.align_msas import *
from compotts.manage_positions import *

class ComPotts_Object:

    def __init__(self, mrf=None, potts_model_file=None, name=None, sequence_file=None, aln_fasta=None, a3m_file=None, input_folder=None, nb_sequences=1000, use_less_sequences=True, hhfilter_threshold=80, perform_filter=True, trimal_gt=0.8, trimal_cons=60, pc_count=1000, reg_lambda_pair_factor=30, trim_alignment=True, rescaling_function="identity", use_w=True, mrf_type=None, hhblits_database=None, min_sequences=1, retry_hhblits_with_memory_limit_if_fail=False, **kwargs):

        self.folder = input_folder

        # MRF
        self.potts_model_file = potts_model_file
        if (self.potts_model_file is None) and (input_folder is not None):
            self.potts_model_file = fm.get_potts_model_file_from_folder(input_folder)
        if mrf is not None:
            self.mrf = mrf
        else:
            if self.potts_model_file is not None:
                self.mrf = Potts_Model.from_msgpack(self.potts_model_file, **kwargs)

        # SEQ_FILE
        if sequence_file is not None:
            self.sequence_file = sequence_file
        elif input_folder is not None:
            self.sequence_file = fm.get_sequence_file_from_folder(input_folder)
        else:
            self.sequence_file = None


        # EXISTING ALIGNMENT FASTA FORMAT
        self.aln_fasta = aln_fasta


        # EXISTING A3M FILE
        self.a3m_file = a3m_file

        # EXISTING FILES
        if input_folder is not None:
            if self.a3m_file is None:
                self.a3m_file = fm.get_a3m_file_from_folder(input_folder)
            self.a3m_reformat = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_reformat.fasta")
            if perform_filter:
                self.filtered = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_filtered_"+str(hhfilter_threshold)+".fasta")
            if use_less_sequences:
                self.less = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_less.fasta")
            if trim_alignment:
                self.trimmed_aln = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_trim_"+str(int(trimal_gt*100))+".fasta")
            if self.aln_fasta is None:
                if self.trimmed_aln is not None:
                    self.aln_fasta = self.trimmed_aln
                elif self.less is not None:
                    self.aln_fasta = self.less
                elif self.filtered is not None:
                    self.aln_fasta = self.filtered
                elif self.a3m_reformat is not None:
                    self.aln_fasta = self.a3m_reformat


        # EXISTING ALIGNMENT FASTA FORMAT
        self.aln_fasta = aln_fasta
#        if self.aln_fasta is None:
#            self.aln_fasta = fm.get_existing_training_set(input_folder, trimal_gt)

        # NAME
        if name is not None:
            self.name = name
        elif input_folder is not None:
            self.name = input_folder.stem
        elif self.sequence_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.sequence_file)
        elif self.aln_fasta is not None:
            self.name = fm.get_first_sequence_clean_name(self.aln_fasta)
        elif self.a3m_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.a3m_file)
        elif mrf is not None:
            self.name = mrf.name
        else:
            self.name = "Billy_"+time.strftime("%Y%m%d-%H%M%S")


        # SEQUENCE
        if self.sequence_file is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.sequence_file)
        elif self.aln_fasta is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.aln_fasta)
        else:
            self.sequence = None



        # WAYS OF TRAINING A POTTS MODEL
        if mrf_type is not None:
            self.mrf_type=mrf_type
        elif self.aln_fasta is not None: # if aln_fasta exists, MRF is trained in a standard way
            self.mrf_type="standard"
        elif self.sequence_file is not None: # if aln_fasta doesn't exist but we have a sequence file, it is used to train the MRF
           self.mrf_type="one_submat"
        else:
            self.mrf_type=mrf_type


        # IF WE NEED TO INFER A POTTS MODEL IN A STANDARD WAY (FROM AN ALIGNMENT)

        if (self.potts_model_file is None) and (self.mrf_type=="standard"):

            # IF WE DON'T HAVE AN A3M FILE AND WE WANT ONE
            if (self.a3m_file is None) and (self.aln_fasta is None):
                self.a3m_file = self.get_folder()/(self.name+".a3m")
                if hhblits_database is None:
                    raise Exception("Must specify a database for hhblits (option -d)")
                else:
                    call_hhblits(self.sequence_file, self.a3m_file, hhblits_database, retry_hhblits_with_memory_limit_if_fail=False, **kwargs)

 
            # REFORMAT A3M_FILE
            if (self.a3m_file is not None) and (self.a3m_reformat is None) and (self.aln_fasta is None):
                self.a3m_reformat = self.get_folder()/(self.name+"_reformat.fasta")
                call_reformat(self.a3m_file, self.a3m_reformat)
                self.aln_fasta = self.a3m_reformat


            # FILTER
            if (self.aln_fasta is not None) and (perform_filter) and (self.filtered is None):
                old_aln_fasta = self.aln_fasta
                self.filtered = self.get_folder()/(self.name+"_filtered_"+str(hhfilter_threshold)+".fasta")
                call_hhfilter(old_aln_fasta, self.filtered, hhfilter_threshold)
                self.aln_fasta = self.filtered


            # USE LESS SEQUENCES
            if (self.aln_fasta is not None) and (use_less_sequences) and (self.less is None):
                old_aln_fasta = self.aln_fasta
                self.less = self.get_folder()/(self.name+"_less.fasta")
                fm.create_fasta_file_with_less_sequences(old_aln_fasta, self.less, nb_sequences)
                self.aln_fasta = self.less


            # TRIM ALIGNMENT
            if (self.aln_fasta is not None) and (trim_alignment) and (self.trimmed_aln is None):
                colnumbering_file = self.get_folder()/(self.name+"_colnumbering.csv")
                old_aln_fasta = self.aln_fasta
                self.trimmed_aln = self.get_folder()/(self.name+"_trim_"+str(int(trimal_gt*100))+".fasta")
                call_trimal(old_aln_fasta, self.trimmed_aln, trimal_gt, trimal_cons, colnumbering_file)
                self.aln_fasta = self.trimmed_aln
                self.real_aln_pos = fm.get_trimal_ncol(colnumbering_file)


            # SET TRAINING SET
            self.training_set = self.aln_fasta
            if self.training_set is None:
                raise Exception('Training alignment missing !')
            if (fm.get_nb_sequences_in_fasta_file(self.training_set)<min_sequences):
                raise Exception("Training set has less than "+str(min_sequences)+" sequences")


            # TRAIN MRF
            self.potts_model_file = (self.get_folder())/(self.name+"_"+self.mrf_type+".mrf")
            self.mrf = Potts_Model.from_training_set(self.training_set, self.potts_model_file, pc_count=pc_count, reg_lambda_pair_factor=reg_lambda_pair_factor, **kwargs)
            if self.mrf is None:
                raise Exception('MRF was not inferred')

            os.system("cp "+str(self.training_set)+" "+str(self.get_folder()/(self.name+"_training_set.fasta")))



        elif (self.potts_model_file is None) and ((self.mrf_type=="one_hot") or (self.mrf_type=="one_submat")):

            # TRAINING SET (SEQUENCE)
            if self.sequence_file is not None:
                self.training_set = self.sequence_file
            elif self.sequence is not None:
                self.sequence_file = self.folder/(self.name+".fasta")
                fm.create_seq_fasta(self.sequence, self.sequence_file, seq_name=self.name)
                self.training_set = self.sequence_file
            else:
                raise Exception('Sequence file missing !')
            self.potts_model_file = (self.get_folder())/(self.name+"_"+self.mrf_type+".mrf")

            # TRAIN MRF
            if (self.mrf_type=="one_hot"):
                self.mrf = Potts_Model.from_sequence_file_to_one_hot(self.training_set, **kwargs)
            elif (self.mrf_type=="one_submat"):
                self.mrf = Potts_Model.from_sequence_file_with_submat(self.training_set, **kwargs)

            os.system("cp "+str(self.training_set)+" "+str(self.get_folder()/(self.name+"_training_set.fasta")))


        else:
            self.training_set = self.aln_fasta 



        # COLNUMBERING
        if trim_alignment:
            colnumbering_file = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_colnumbering.csv")
            self.real_aln_pos = fm.get_trimal_ncol(colnumbering_file)
        elif self.aln_fasta is not None:
            nb_pos = fm.get_nb_columns_in_alignment(self.aln_fasta) 
            self.real_aln_pos = [pos for pos in range(nb_pos)]



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


        # RESCALING
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

def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--potts_model_file', help="Potts model", type=pathlib.Path)
    parser.add_argument('-s', '--sequence_file', help="Sequence file", type=pathlib.Path)
    parser.add_argument('-a3m', '--a3m_file', help="HH-blits output file", type=pathlib.Path)
    parser.add_argument('-f', '--input_folder', help="Folder containing files for sequence", type=pathlib.Path)
    parser.add_argument('-aln', '--aln_fasta', help="Alignment file in fasta format", type=pathlib.Path)
    parser.add_argument('-n', '--nb_sequences', help="Max number of sequences in the MRF training alignment", default=1000, type=int)
    parser.add_argument('-nmin', '--min_sequences', help="Min number of sequences in the MRF training alignment", default=1, type=int)
    parser.add_argument('-m', '--mrf_type', help="Mode", choices=('standard', 'one_hot', 'one_submat'), default='standard')
    parser.add_argument('-nm', '--name', help="Name")

    # trimal
    parser.add_argument('-trimgt', '--trimal_gt', help="trimal gt", default=0.8, type=float)
    parser.add_argument('-trimcons', '--trimal_cons', help="trimal cons", default=60, type=float)


    # CCMpredPy options
    parser.add_argument('--pc_count', help="CCMpred : Specify number of pseudocounts (default : 1000)", default=1000)
    parser.add_argument('--reg_lambda_pair_factor', help="CCMpred : Regularization parameter for pair potentials (L2 regularization with lambda_pair = lambda_pair-factor * scaling) [CCMpred default: 0.2, our default : 30]", default=30)

    args = vars(parser.parse_args(args))

    obj = ComPotts_Object(**args)


if __name__=="__main__":
    main()
