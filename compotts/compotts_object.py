""" ComPotts Object : contains MSA files locations, Potts models, etc. """

import os
import time
import argparse
import sys
import shutil
import basic_modules.files_management as fm
from basic_modules.util import *
from basic_modules.tool_wrapper import *
from basic_modules.potts_model import *
from compotts.rescaling import *
from compotts.align_msas import *
from compotts.manage_positions import *
from compotts.find_cutoff_index import *

class ComPotts_Object:

    def __init__(self, mrf=None, potts_model_file=None, input_folder=None, mrf_type=None, sequence_file=None, aln_fasta=None, sequences_fetcher='hhblits', a3m_file=None, name=None, hhblits_database=None, blast_fasta=None, use_evalue_cutoff = False, hhr_file=None, blast_xml=None, filter_alignment=True, hhfilter_threshold=80, use_less_sequences=True, max_nb_sequences=1000, min_nb_sequences=1, trim_alignment=True, trimal_gt=0.8, trimal_cons=0, pc_count=None, reg_lambda_pair_factor=30, rescaling_function="identity", **kwargs):

        self.folder = input_folder
        self.mrf = mrf
        self.colnumbering_file = None

        # POTTS MODEL
        self.potts_model_file = potts_model_file
        if (self.potts_model_file is None) and (input_folder is not None):
            self.potts_model_file = fm.get_potts_model_file_from_folder(input_folder, mrf_type=mrf_type)
        if mrf is not None:
            self.mrf = mrf
        else:
            if self.potts_model_file is not None:
                if os.path.isfile(self.potts_model_file):
                    self.mrf = Potts_Model.from_msgpack(self.potts_model_file, **kwargs)

        # SEQUENCE FILE
        if sequence_file is not None:
            self.sequence_file = sequence_file
        elif input_folder is not None:
            self.sequence_file = fm.get_sequence_file_from_folder(input_folder)
        else:
            self.sequence_file = None
        

        # EXISTING TRAINING SET ?
        self.training_set = None
        if (self.folder is not None):
            self.training_set = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_training_set.fasta")


        # NAME
        if name is not None:
            self.name = name
        elif input_folder is not None:
            self.name = input_folder.stem
        elif self.sequence_file is not None:
            self.name = fm.get_first_sequence_clean_name(self.sequence_file)
        elif aln_fasta is not None:
            self.name = fm.get_first_sequence_clean_name(aln_fasta)
        elif self.mrf is not None:
            self.name = self.mrf.name
        elif a3m_file is not None:
            self.name = fm.get_first_sequence_clean_name(a3m_file)
        elif self.training_set is not None:
            self.name = fm.get_first_sequence_clean_name(self.training_set)
        else:
            self.name = "Billy_"+time.strftime("%Y%m%d-%H%M%S")


        # CHECK IF EXISTING A3M OR BLAST FASTA
        if a3m_file is not None:
            self.a3m_file = a3m_file
        elif self.folder is not None:
            self.a3m_file = fm.get_a3m_file_from_folder(self.folder)

        self.reformat = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_reformat.fasta")

        if blast_fasta is not None:
            self.blast_fasta = blast_fasta
        elif self.folder is not None:
            self.blast_fasta = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_blast.fasta")


        # WAYS OF TRAINING A MRF
        if mrf_type is not None:
            self.mrf_type = mrf_type
        elif (aln_fasta is not None) or (self.a3m_file is not None) or (self.blast_fasta is not None):
            self.mrf_type = "standard"
        elif self.sequence_file is not None:
            self.mrf_type = "one_submat"
        else: # default
            self.mrf_type = "standard"


        # TRAINING SET DOES NOT EXIST
        if (self.training_set is None):

            if aln_fasta is not None:
                self.training_set = aln_fasta

            elif (mrf is None): # AND WE NEED ONE

                if (self.mrf_type=="standard"):

                    if sequences_fetcher=='hhblits': # IF WE USE HHBLITS TO GET THE ALIGNMENT
                        if self.a3m_file is None:
                            self.a3m_file = self.get_folder()/(self.name+".a3m")
                            if hhblits_database is not None:
                                call_hhblits(self.sequence_file, self.a3m_file, hhblits_database, retry_hhblits_with_memory_limit_if_fail=False, **kwargs)
                            else:
                                raise Exception("No HHblits database (-d option)")
                        self.reformat = self.get_folder()/(self.name+"_reformat.fasta")
                        call_reformat(self.a3m_file, self.reformat)
                        self.training_set = self.reformat


                    elif sequences_fetcher=='blast':
                        if self.blast_fasta is None:
                            raise Exception("Need BLAST fasta output (BLAST call not implemented yet)")
                        self.training_set = self.blast_fasta


                    if use_evalue_cutoff:
                        if sequences_fetcher=='hhblits':
                            if hhr_file is not None:
                                self.hhr_file = hhr_file
                            else:
                                self.hhr_file = fm.get_file_from_folder_ending_with_extension(self.folder(), ".hhr")
                            if self.hhr_file is None:
                                raise Exception("Need a .hhr file !")
                            else:
                                cutoff_index = find_hhblits_cutoff_index(self.hhr_file)
                        elif sequences_fetcher=='blast':
                            if blast_xml is not None:
                                self.blast_xml = blast_xml
                            else:
                                self.blast_xml = fm.get_file_from_folder_ending_with_extension(self.folder(), "_blast.xml")
                            if self.blast_xml is None:
                                raise Exception("Need a BLAST XML output file !")
                            else:
                                cutoff_index = find_blast_cutoff_index(self.blast_xml)
                        else:
                            raise Exception("Unknown sequences fetcher")

                        self.cutoff_fasta = self.get_folder()/(self.name+"_cutoff_"+str(cutoff_index)+".fasta")
                        fm.create_fasta_file_with_less_sequences(self.training_set, self.cutoff_fasta, cutoff_index)
                        self.training_set = self.cutoff_fasta


                    if filter_alignment:
                        self.filtered = self.get_folder()/(self.name+"_filtered_"+str(hhfilter_threshold)+".fasta")
                        call_hhfilter(self.training_set, self.filtered, hhfilter_threshold)
                        self.training_set = self.filtered

                    if use_less_sequences:
                        self.less = self.get_folder()/(self.name+"_less.fasta")
                        fm.create_fasta_file_with_less_sequences(self.training_set, self.less, max_nb_sequences)
                        self.training_set = self.less

                    if trim_alignment:
                        self.colnumbering_file = self.get_folder()/(self.name+"_colnumbering.csv")
                        self.trimmed_aln = self.get_folder()/(self.name+"_trim_"+str(int(trimal_gt*100))+".fasta")
                        call_trimal(self.training_set, self.trimmed_aln, trimal_gt, trimal_cons, self.colnumbering_file)
                        self.training_set = self.trimmed_aln

                elif (self.mrf_type=="one_submat") or (self.mrf_type=="one_hot"):
                    if self.sequence_file is None:
                        raise Exception("Need a sequence file for one_submat training")
                    else:
                        self.training_set = self.sequence_file
                else:
                    raise Exception("Unknown MRF training type")

            shutil.copy(str(self.training_set), str(self.get_folder()/(self.name+"_training_set.fasta")))


        if (self.mrf is None): # WE NEED TO TRAIN AN MRF

            if potts_model_file is not None:
                self.potts_model_file = potts_model_file
            else:
                self.potts_model_file = (self.get_folder())/(self.name+"_"+self.mrf_type+".mrf")
            

            if fm.get_nb_sequences_in_fasta_file(self.training_set)<min_nb_sequences:
                raise Exception("Nb sequences in the training set is less than "+str(nb_sequences))

            if self.mrf_type=="standard":
                if pc_count is None:
                    pc_count = fm.get_nb_sequences_in_fasta_file(self.training_set) # DEFAULT : NB PSEUDOCOUNTS = NB SEQUENCES
                self.mrf = Potts_Model.from_training_set(self.training_set, self.potts_model_file, pc_count=pc_count, reg_lambda_pair_factor=reg_lambda_pair_factor, **kwargs)

            elif self.mrf_type=="one_hot":
                self.mrf = Potts_Model.from_sequence_file_to_one_hot(self.training_set, filename=self.potts_model_file, **kwargs)

            elif self.mrf_type=="one_submat":
                self.mrf = Potts_Model.from_sequence_file_with_submat(self.training_set, filename=self.potts_model_file, **kwargs)


            self.original_mrf = self.mrf
            if (rescaling_function!="identity"):
                print("rescaling Potts model")
                self.mrf = get_rescaled_mrf(self.mrf, rescaling_function, use_w=use_w)


        if self.mrf is None:
            raise Exception("Potts model could not be inferred.")


        # COLNUMBERING : TRAINING SET POSITIONS -> "REFORMAT"/"BLAST FASTA" POSITIONS
        if trim_alignment:
            if self.colnumbering_file is None:
                self.colnumbering_file = fm.get_file_from_folder_ending_with_extension(self.get_folder(), "_colnumbering.csv")
            if self.colnumbering_file is not None:
                self.real_aln_pos = fm.get_trimal_ncol(self.colnumbering_file)
            else:
                nb_pos = fm.get_nb_columns_in_alignment(self.training_set)
                self.real_aln_pos = [pos for pos in range(nb_pos)]
        else:
            nb_pos = fm.get_nb_columns_in_alignment(self.training_set)
            self.real_aln_pos = [pos for pos in range(nb_pos)]

        # SEQUENCE
        if self.sequence_file is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(self.sequence_file)
        elif aln_fasta is not None:
            self.sequence = fm.get_first_sequence_in_fasta_file(aln_fasta)
        elif (self.a3m_file is not None) and (sequences_fetcher=='hhblits'):
            self.sequence = fm.get_first_sequence_in_fasta_file(self.a3m_file)
        elif (self.blast_fasta is not None) and (sequences_fetcher=='blast'):
            self.sequence = fm.get_first_sequence_in_fasta_file(self.blast_fasta)
        elif (self.training_set is not None):
            self.sequence = fm.get_first_sequence_in_fasta_file(self.training_set)
        else:
            self.sequence = None

        # "REFORMAT"/"BLAST FASTA" POSITIONS -> SEQUENCE POSITIONS
        if (self.sequence is not None) and (self.training_set is not None):
            if (sequences_fetcher=='blast') and (self.blast_fasta is not None):
                aln_first_seq = fm.get_first_sequence_in_fasta_file(self.blast_fasta)
                seq_aln_pos = get_real_pos_list(self.sequence, aln_first_seq)
                self.real_seq_pos = [seq_aln_pos[pos] for pos in self.real_aln_pos]
            elif (sequences_fetcher=='hhblits') and (self.reformat is not None):
                aln_first_seq = fm.get_first_sequence_in_fasta_file(self.reformat)
                seq_aln_pos = get_real_pos_list(self.sequence, aln_first_seq)
                self.real_seq_pos = [seq_aln_pos[pos] for pos in self.real_aln_pos]
            else:
                self.real_seq_pos = self.real_aln_pos


    def get_folder(self):
        if not hasattr(self, "folder"):
            self.folder = fm.create_foler(self.name)
        elif self.folder is None:
            self.folder = fm.create_folder(self.name)
        if (not self.folder.is_dir()):
            self.folder = fm.create_folder(self.folder)
        return self.folder

    def get_seq_positions(self, positions):
        return [self.real_seq_pos[pos] for pos in positions]

    #@classmethod
    #def from_merge(cls, obj1, obj2, aligned_positions, output_folder, **kwargs):
    #TODO

def main(args=sys.argv[1:]):

    parser = argparse.ArgumentParser()
    parser.add_argument('-nm', '--name', help="Name")

    # existing input
    parser.add_argument('-p', '--potts_model_file', help="Potts model", type=pathlib.Path)
    parser.add_argument('-s', '--sequence_file', help="Sequence file", type=pathlib.Path)
    parser.add_argument('-f', '--input_folder', help="Folder containing files for sequence", type=pathlib.Path)
    parser.add_argument('-m', '--mrf_type', help="Mode", choices=('standard', 'one_hot', 'one_submat'), default='standard')
    parser.add_argument('-aln', '--aln_fasta', help="Alignment file in fasta format", type=pathlib.Path)
    parser.add_argument('-a3m', '--a3m_file', help="HH-blits .a3m output file", type=pathlib.Path)
    parser.add_argument('-blastf', '--blast_fasta', help="BLAST fasta alignment output file", type=pathlib.Path) # TODO BLAST FASTA -> MUSCLE

    # fetch sequences
    parser.add_argument('--sequences_fetcher', help="How we get our sequences (hhblits or blast)", default="hhblits")
    parser.add_argument('-d', '--hhblits_database', help="HHblits database path", type=pathlib.Path)

    # use E-value cutoff
    parser.add_argument('-cutev', '--use_evalue_cutoff', help="Use E-value cutoff", default=False, action='store_true')
    parser.add_argument('-hhr', '--hhr_file', help="HH-blits .hhr output file", type=pathlib.Path)
    parser.add_argument('-blastxml', '--blast_xml', help="BLAST XML output file", type=pathlib.Path)


    # hhfilter options
    parser.add_argument('-nofilt', '--dont_filter_alignment', help="don't filter alignment with hhfilter (default = do)", default=False, action='store_true')
    parser.add_argument('-hht', '--hhfilter_threshold', help="hhfilter_threshold", default=80, type=float)
    
    # use less sequences
    parser.add_argument('--dont_use_less_sequences', help="Don't limit the number of sequences in the training set (default = limit to 1000)", default=False, action='store_true')
    parser.add_argument('-nmax', '--max_nb_sequences', help="Max number of sequences in the MRF training alignment", default=1000, type=int)
    parser.add_argument('-nmin', '--min_nb_sequences', help="Min number of sequences in the MRF training alignment (raises an Exception if nb sequences in the training set < nmin)", default=1, type=int)

    # trimal
    parser.add_argument('-notrim', '--dont_trim_alignment', help="don't perform trim", default=False, action='store_true') 
    parser.add_argument('-trimgt', '--trimal_gt', help="trimal gt", default=0.8, type=float)
    parser.add_argument('-trimcons', '--trimal_cons', help="trimal cons", default=0, type=float)


    # CCMpredPy options
    parser.add_argument('--pc_count', help="CCMpred : Specify number of pseudocounts (default : 1000)", default=1000)
    parser.add_argument('--reg_lambda_pair_factor', help="CCMpred : Regularization parameter for pair potentials (L2 regularization with lambda_pair = lambda_pair-factor * scaling) [CCMpred default: 0.2, our default : 30]", default=30)

    # rescaling
    parser.add_argument('--rescaling_function', help="Rescaling function for the inferred Potts model (default : none (identity))", default="identity")

    args = vars(parser.parse_args(args))
    args["use_less_sequences"] = not args["dont_use_less_sequences"]
    args["trim_alignment"] = not args["dont_trim_alignment"]
    args["filter_alignment"] = not args["dont_filter_alignment"]

    obj = ComPotts_Object(**args)


if __name__=="__main__":
    main()
