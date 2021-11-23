""" Potts Objects are objects handled by PPalign, including Potts models, sequences, MSAs and correspondences between them """
import uuid
import sys
import argparse
import copy

from comutils.tool_wrapper import *
from comutils.find_cutoff_index import *
from comutils.blast_utils import *

from ppalign.manage_positions import *
from ppalign.compute_scores import *

from makepotts.rescaling import *
from makepotts.potts_model import *
#from makepotts.handle_insertions import *
from infer_insertion_penalties.__main__ import *

class Potts_Object:

    def __init__(self):
        pass


    @classmethod
    def from_folder(cls, potts_folder, v_rescaling_function="identity", w_rescaling_function="identity", use_w=True, use_insertion_penalties=False, **kwargs):
        """ instantiate Potts object from Potts folder """
        feature = cls()

        feature.folder = potts_folder
        if not feature.folder.is_dir():
            raise Exception("Feature folder does not exist")

        if (potts_folder/"aln_train.fasta").is_file():
            feature.aln_train = potts_folder/"aln_train.fasta"
        else:
            feature.aln_train = None

        if (potts_folder/"aln_before_trim.fasta").is_file():
            feature.aln_before_trim = potts_folder/"aln_before_trim.fasta"
        else:
            feature.aln_before_trim = None
        

        if (potts_folder/"sequence.fasta").is_file():
            feature.sequence_file = potts_folder/"sequence.fasta"
            feature.sequence = fm.get_first_sequence_in_fasta_file(feature.sequence_file)
        else:
            feature.sequence = None

        feature.potts_model = None
        try:
            feature.potts_model_file = potts_folder/"potts_model.mrf"
            feature.potts_model = Potts_Model.from_msgpack(feature.potts_model_file)
        except Exception as e:
            print("Potts model was not found")
            feature.potts_model = None

        feature.mrf_pos_to_seq_pos=None
        try:
            feature.mrf_pos_to_seq_pos = fm.get_list_from_csv(potts_folder/"mrf_pos_to_seq_pos.csv") # mrf_pos_to_aln_pos[i] = position in sequence corresponding to position i in Potts model
        except Exception as e:
            feature.mrf_pos_to_seq_pos = None

        feature.aln_pos_to_seq_pos=None
        try:
            feature.aln_pos_to_seq_pos = fm.get_list_from_csv(potts_folder/"aln_pos_to_seq_pos.csv") # aln_pos_to_aln_pos[i] = position in sequence corresponding to position i in original MSA
        except Exception as e:
            feature.aln_pos_to_seq_pos = None

        feature.mrf_pos_to_aln_pos=None
        try:
            feature.mrf_pos_to_aln_pos = fm.get_list_from_csv(potts_folder/"mrf_pos_to_aln_pos.csv") # mrf_pos_to_aln_pos[i] = position in original_aln corresponding to position i in Potts model
        except Exception as e:
            feature.mrf_pos_to_aln_pos = None

        if (feature.potts_model is not None):
            feature.potts_model = get_rescaled_potts_model(feature.potts_model, v_rescaling_function, w_rescaling_function, use_w=use_w, **kwargs)

        if (potts_folder/"insertion_penalties.tsv").is_file() and use_insertion_penalties:
            feature.insertion_penalties = get_insertion_penalties_from_file(potts_folder/"insertion_penalties.tsv")
            assert(len(feature.insertion_penalties['open']) == len(feature.insertion_penalties['extend']))
            assert(len(feature.insertion_penalties['open']) == feature.potts_model.ncol+1)
        else:
            feature.insertion_penalties = None
        
        return feature


    


    @classmethod
    def from_sequence_alone(cls, potts_folder, sequence_file, inference_type, **kwargs):

        if not potts_folder.is_dir():
            potts_folder.mkdir()

        potts_model_file = potts_folder/"potts_model.mrf"
        if inference_type=="one_submat":
            potts_model = Potts_Model.from_sequence_file_with_submat(sequence_file, filename=potts_model_file, **kwargs)
        elif inference_type=="one_hot":
            potts_model = Potts_Model.from_sequence_file_to_one_hot(sequence_file, filename=potts_model_file, **kwargs)
        else:
            raise Exception("Unknown inference type")

        nb_pos = fm.get_nb_columns_in_alignment(sequence_file)
        mrf_pos_to_aln_pos = [pos for pos in range(nb_pos)]
        mrf_pos_to_seq_pos = [pos for pos in range(nb_pos)]
        aln_pos_to_seq_pos = [pos for pos in range(nb_pos)]

        return cls.from_potts_model(potts_folder, potts_model_file, sequence_file=sequence_file, aln_train=sequence_file, aln_before_trim=sequence_file, mrf_pos_to_aln_pos=mrf_pos_to_aln_pos, aln_pos_to_seq_pos=aln_pos_to_seq_pos, mrf_pos_to_seq_pos=mrf_pos_to_seq_pos, **kwargs)


    @classmethod
    def from_sequence_with_blast(cls, potts_folder, sequence_file, blast_database, nb_sequences_blast=100000, blast_evalue=1, **kwargs):
        unaligned_fasta = potts_folder/"unaligned_sequences.fasta"
        blast_xml = potts_folder/"blast.xml"
        blast_xml, unaligned_fasta = get_blast_xml_and_fasta_output_from_sequence_file(sequence_file, blast_database, blast_fasta=unaligned_fasta, blast_xml=blast_xml, n=nb_sequences_blast, evalue=blast_evalue)
        return cls.from_blast_files(potts_folder, sequence_file, unaligned_fasta, blast_xml, **kwargs)


    @classmethod
    def from_blast_files(cls, potts_folder, sequence_file, unaligned_fasta, blast_xml=None, use_evalue_cutoff=False, use_less_sequences=True, max_nb_sequences=1000, **kwargs):


        if (use_evalue_cutoff) or (use_less_sequences):
            if use_evalue_cutoff:
                if blast_xml is None:
                    raise Exception("Need BLAST XML file to determine cutoff index")
                max_nb_sequences = find_blast_cutoff_index(blast_xml)
            less_file = potts_folder/("less.a3m")
            fm.create_file_with_less_sequences(unaligned_fasta, less_file, max_nb_sequences)
            unaligned_fasta = less_file

        fm.add_sequence_to_fasta_file_if_missing(unaligned_fasta, sequence_file)
        clean_unaligned_fasta = potts_folder/"clean_unaligned_sequences.fasta"
        fm.remove_sequences_with_bad_characters(unaligned_fasta, clean_unaligned_fasta)
        tmp_unaligned_fasta = potts_folder/"clean_unaligned_sequences_bis.fasta"
        fm.copy(clean_unaligned_fasta, tmp_unaligned_fasta)
        fm.create_file_with_less_sequences(tmp_unaligned_fasta, clean_unaligned_fasta, nb_sequences=10000, fileformat="fasta") # limit to 10000 sequences for MAFFT
        tmp_unaligned_fasta.unlink()
        aln_mafft = potts_folder/"aln_mafft.fasta"
        call_mafft(clean_unaligned_fasta, aln_mafft)
        aln_file = potts_folder/"aln_clean.fasta"
        fm.remove_positions_with_gaps_in_first_sequence(aln_mafft, aln_file) #TODO -> insertions
        return cls.from_aln_file(potts_folder, aln_file, sequence_file=sequence_file, **kwargs)



    @classmethod
    def from_sequence_with_hhblits(cls, potts_folder, sequence_file, hhblits_database, **kwargs):
        aln_with_insertions = potts_folder/"hhblits_output.a3m"
        hhr_file = call_hhblits(sequence_file, aln_with_insertions, hhblits_database, **kwargs)
        return cls.from_hhblits_files(potts_folder, aln_with_insertions, sequence_file=sequence_file, hhr_file=hhr_file, **kwargs)




    @classmethod
    def from_hhblits_files(cls, potts_folder, aln_with_insertions, sequence_file=None, hhr_file=None, use_evalue_cutoff=False, use_less_sequences=True, max_nb_sequences=1000, filter_alignment=True, hhfilter_threshold=80, **kwargs):

        clean_file = potts_folder/"clean.a3m"
        fm.remove_sequences_with_bad_characters(aln_with_insertions, clean_file)
        aln_with_insertions = clean_file

        if filter_alignment:
            filtered = potts_folder/"filtered.a3m"
            call_hhfilter(aln_with_insertions, filtered, hhfilter_threshold)
            aln_with_insertions = filtered
 
        if (use_evalue_cutoff) or (use_less_sequences):
            if use_evalue_cutoff:
                if hhr_file is None:
                    raise Exception("Need hhr file to determine cutoff index")
                max_nb_sequences = find_hhblits_cutoff_index(hhr_file)
            less_file = potts_folder/("less.a3m")
            fm.create_file_with_less_sequences(aln_with_insertions, less_file, max_nb_sequences)
            aln_with_insertions = less_file

        aln_file = potts_folder/"reformat.fasta"
        call_reformat(aln_with_insertions, aln_file)



        # CREATE SEQUENCE FILE IF NONE
        if sequence_file is None:
            seq = fm.get_first_sequence_in_fasta_file(aln_with_insertions)
            seq_name = fm.get_first_sequence_name(aln_with_insertions)
            sequence_file = potts_folder/"sequence.fasta"
            fm.create_seq_fasta(seq, sequence_file, seq_name=seq_name)

        return cls.from_aln_file(potts_folder, aln_file, sequence_file=sequence_file, aln_with_insertions=aln_with_insertions, **kwargs)





    @classmethod
    def from_aln_file(cls, potts_folder, aln_file, sequence_file=None, trim_alignment=False, trimal_gt=0.8, trimal_cons=0, min_nb_sequences=1, max_potts_model_length=250, pc_single_count=1, aln_with_insertions=None, infer_potts_model=True, potts_model_file=None, **kwargs):

        aln_train = aln_file

        # TRIM IF NEEDED
        if trim_alignment:
            colnumbering_file = potts_folder/"colnumbering.csv"
            trimmed_aln = potts_folder/("trim.fasta")
            call_trimal(aln_train, trimmed_aln, trimal_gt, trimal_cons, colnumbering_file)
            aln_train = trimmed_aln
            mrf_pos_to_aln_pos = fm.get_trimal_ncol(colnumbering_file)
        else:
            nb_pos = fm.get_nb_columns_in_alignment(aln_train)
            mrf_pos_to_aln_pos = [pos for pos in range(nb_pos)]


        # CHECK ALN TRAIN OK
        if fm.get_nb_sequences_in_fasta_file(aln_train)<min_nb_sequences:
                raise Exception("Less than "+str(min_nb_sequences)+" in the training set : "+str(fm.get_nb_sequences_in_fasta_file(aln_train)))

        if fm.get_nb_columns_in_alignment(aln_train)>max_potts_model_length:
            raise Exception("More than "+str(max_potts_model_length)+" columns in the alignment, won't infer the Potts model.")

        # INFER POTTS MODEL
        if infer_potts_model:
            potts_model_file = potts_folder/"potts_model.mrf"
            potts_model = Potts_Model.from_training_set(aln_train, potts_model_file, pc_single_count=pc_single_count, **kwargs)
        elif potts_model_file is not None:
            # if Potts model is not inferred and from elsewhere, copy file
            new_potts_model_file = potts_folder/"potts_model.mrf"
            if not new_potts_model_file.is_file():
                fm.copy(potts_model_file, new_potts_model_file)
            potts_model_file = new_potts_model_file


       # HANDLE SEQUENCE FILE
        if sequence_file is not None:
            original_first_seq = fm.get_first_sequence_in_fasta_file(aln_file)
            seq = fm.get_first_sequence_in_fasta_file(sequence_file)
            mrf_pos_to_seq_pos = get_mrf_pos_to_seq_pos(original_first_seq, seq, mrf_pos_to_aln_pos)
            aln_pos_to_seq_pos = get_pos_first_seq_to_second_seq(original_first_seq, seq) 
        else: # if no sequence file is provided, sequence is the first sequence of the MSA
            seq = fm.get_first_sequence_in_fasta_file(aln_file)
            seq_name = fm.get_first_sequence_name(aln_file)
            sequence_file = potts_folder/"sequence.fasta"
            fm.create_seq_fasta(seq, sequence_file, seq_name=seq_name)
            mrf_pos_to_seq_pos = mrf_pos_to_aln_pos
            aln_pos_to_seq_pos = [pos for pos in range(len(seq))]

        

        return cls.from_potts_model(potts_folder, potts_model_file, sequence_file=sequence_file, aln_train=aln_train, aln_before_trim=aln_file, aln_with_insertions=aln_with_insertions, aln_pos_to_seq_pos=aln_pos_to_seq_pos, mrf_pos_to_seq_pos=mrf_pos_to_seq_pos, mrf_pos_to_aln_pos=mrf_pos_to_aln_pos, **kwargs)




    @classmethod
    def from_potts_model(cls, potts_folder, potts_model_file, v_rescaling_function="identity", w_rescaling_function="identity", sequence_file=None, aln_train=None, aln_before_trim=None, aln_with_insertions=None, mrf_pos_to_aln_pos=None, aln_pos_to_seq_pos=None, mrf_pos_to_seq_pos=None, insert_null_at_trimmed=False, insert_v_star_at_trimmed=False, v_null_is_v0=True, use_insertion_penalties=False, keep_tmp_files=False, pc_insertions_tau=0, light=False, **kwargs):

        if potts_model_file is not None:
            potts_model = Potts_Model.from_msgpack(potts_model_file)


            if ((insert_null_at_trimmed) or (insert_v_star_at_trimmed)): # RE-INSERT NULL COLUMNS
                if aln_before_trim is None:
                    raise Exception("MSA before trim should be provided")
                if mrf_pos_to_aln_pos is None:
                    raise Exception("mrf_pos_to_aln_pos is None")
                if (insert_v_star_at_trimmed):
                    potts_model.insert_vi_star_gapped_to_complete_mrf_pos(mrf_pos_to_aln_pos, fm.get_nb_columns_in_alignment(aln_before_trim), aln_before_trim)
                elif (insert_null_at_trimmed):
                    if v_null_is_v0:
                        v_null = np.tile(get_background_v0("identity", rescale_removed_v0=False), (1,1))
                    else:
                        v_null = np.zeros((1,21))
                    potts_model.insert_null_positions_to_complete_mrf_pos(mrf_pos_to_aln_pos, fm.get_nb_columns_in_alignment(aln_before_trim), v_null=v_null)
                mrf_pos_to_seq_pos = aln_pos_to_seq_pos
                mrf_pos_to_aln_pos = [pos for pos in range(fm.get_nb_columns_in_alignment(aln_before_trim))]
                potts_model.to_msgpack(potts_model_file)
            elif use_insertion_penalties and ((aln_with_insertions is not None) or (aln_before_trim is not None)) and (mrf_pos_to_aln_pos is not None): # LOWER CASE FOR INSERTION PENALTIES
                if aln_with_insertions is None:
                    aln_with_insertions = aln_before_trim
                aln_with_insertions_and_trim = potts_folder/"aln_with_insertions_and_trim.a3m"
                columns_not_trimmed = mrf_pos_to_aln_pos
                lower_case_trimmed_columns(aln_with_insertions, aln_with_insertions_and_trim, columns_not_trimmed)
                aln_with_insertions = aln_with_insertions_and_trim

            if (v_rescaling_function!="identity") and (w_rescaling_function!="identity"):
                potts_model = get_rescaled_potts_model(potts_model, v_rescaling_function, w_rescaling_function, use_w=use_w, **kwargs)
                potts_model.to_msgpack(potts_model_file)


        # INSERTION PENALTIES
        if use_insertion_penalties:
            if aln_with_insertions is None:
                raise Exception("File with insertions as lower case is required")
            insertions_file = potts_folder/"insertion_penalties.tsv"
            infer_insertion_penalties_in_file(aln_with_insertions, insertions_file, pc_insertions_tau=pc_insertions_tau)


        # HANDLE CORRESPONDENCES BETWEEN POTTS MODEL POSITIONS AND SEQ/MSA POSITIONS
        if mrf_pos_to_seq_pos is not None:
            fm.write_list_to_csv(mrf_pos_to_seq_pos, potts_folder/"mrf_pos_to_seq_pos.csv")
        if mrf_pos_to_aln_pos is not None:
            fm.write_list_to_csv(mrf_pos_to_aln_pos, potts_folder/"mrf_pos_to_aln_pos.csv")
        if aln_pos_to_seq_pos is not None:
            fm.write_list_to_csv(aln_pos_to_seq_pos, potts_folder/"aln_pos_to_seq_pos.csv")


        # COPY FILES IN POTTS FOLDER
        if sequence_file is not None:
            fm.copy(sequence_file, potts_folder/"sequence.fasta")
        if aln_train is not None:
            fm.copy(aln_train, potts_folder/"aln_train.fasta")
        if aln_before_trim is not None:
            fm.copy(aln_before_trim, potts_folder/"aln_before_trim.fasta")
        if aln_with_insertions is not None:
            fm.copy(aln_with_insertions, potts_folder/"aln_with_insertions.fasta")
        if potts_model_file is not None:
            fm.copy(potts_model_file, potts_folder/"potts_model.mrf")

        # REMOVE TEMPORARY FILES
        if not keep_tmp_files:
            for fn in ["hhblits_output.a3m", "hhblits_output.hhr", "filtered.a3m", "less.a3m", "reformat.fasta", "trim.fasta", "colnumbering.csv", "aln_with_insertions_and_trim.a3m", "clean.a3m"]:
                if (potts_folder/fn).is_file():
                    fm.remove_file(potts_folder/fn)

        if light:
            for fn in ["aln_train.fasta", "aln_with_insertions.fasta", "aln_before_trim.fasta", "potts_model_mrf_README.txt"]:
                if (potts_folder/fn).is_file():
                    fm.remove_file(potts_folder/fn)

        return cls.from_folder(potts_folder, v_rescaling_function="identity", w_rescaling_function="identity", use_insertion_penalties=use_insertion_penalties)




    @classmethod
    def from_merge(cls, potts_folder, objects, aligned_positions_dict, use_less_sequences=False, **kwargs):
        """ create Potts object by merging two Potts objects aligned """
        if not potts_folder.is_dir():
            potts_folder.mkdir()
        seed_aln = potts_folder/"seed_aln.fasta"
        get_original_msas_aligned_from_aligned_positions(aligned_positions_dict, objects, seed_aln)
        return cls.from_folder(potts_folder=potts_folder, aln_file=seed_aln, use_less_sequences=use_less_sequences, **kwargs)



    def get_name(self):
        if self.sequence_file is not None:
            self.name = fm.get_first_sequence_name(self.sequence_file)
        elif self.seed_aln is not None:
            self.name = fm.get_first_sequence_name(self.seed_aln)
        elif self.aln_train is not None:
            self.name = fm.get_first_sequence_name(self.aln_train)
        elif self.potts_model is not None:
            self.name = fm.get_first_sequence_name(self.potts_model)
        elif self.potts_model_file is not None:
            potts_model = Potts_Model.from_msgpack(potts_model_file)
            self.name = fm.get_first_sequence_name(self.potts_model)
        else:
            self.name = None
        return self.name


    def get_seq_positions(self, positions):
        """ seq_positions[i] is the position in the sequence corresponding to position i in the Potts model """
        seq_positions = []
        for pos in positions:
            if pos is None:
                seq_positions.append(None)
            elif pos is '-':
                seq_positions.append('-')
            else:
                seq_positions.append(self.mrf_pos_to_seq_pos[pos])
        return seq_positions

    def get_aln_positions(self, positions):
        """ aln_positions[i] is the position in the original MSA corresponding to position i in the Potts model """
        aln_positions = []
        for pos in positions:
            if pos is None:
                aln_positions.append(None)
            else:
                aln_positions.append(self.mrf_pos_to_aln_pos[pos])
        return aln_positions


    def get_seq_pos_to_mrf_pos(self):
        """ seq_pos_to_mrf_pos[i] is the position in the Potts model corresponding to position i in the sequence """
        seq_pos_to_mrf_pos = []
        for pos in range(len(self.sequence)):
            try:
                mrf_pos = self.mrf_pos_to_seq_pos.index(pos)
            except Exception as e:
                mrf_pos = None
            seq_pos_to_mrf_pos.append(mrf_pos)
        return seq_pos_to_mrf_pos


    def get_positions_in_sequence_that_are_not_in_train_alignment(self):
        """ returns all positions that are in the sequence but not in the MSA used to train the Potts model """
        seq_pos_to_mrf_pos = self.get_seq_pos_to_mrf_pos()
        in_seq_not_in_aln = []
        for pos in range(len(self.sequence)):
            if seq_pos_to_mrf_pos is None:
                in_seq_not_in_aln.append(seq_pos_to_mrf_pos)
        return in_seq_not_in_aln


    def insert_null_at_trimmed(self, remove_v0=False, change_mrf_pos_lists=False, **kwargs):
        """ inserts null columns at positions trimmed by trimal """
        if remove_v0:
            v_null = np.tile(get_background_v0(**kwargs), (1,1))
        else:
            v_null = np.zeros((1,21)) 
        self.potts_model.insert_null_positions_to_complete_mrf_pos(self.mrf_pos_to_aln_pos, len(self.sequence), v_null=v_null)
        if change_mrf_pos_lists:
            self.mrf_pos_to_seq_pos = aln_pos_to_seq_pos
            self.mrf_pos_to_aln_pos = [pos for pos in range(fm.get_nb_columns_in_alignment(seed_aln))]


    def to_folder(self, output_folder):
        if not output_folder.is_dir():
            output_folder.mkdir()
        self.potts_model.to_msgpack(output_folder/"potts_model.mrf")
        fm.copy(self.seed_aln, output_folder/"seed_aln.fasta")
        fm.copy(self.aln_train, output_folder/"aln_train.fasta")
        fm.create_seq_fasta(self.sequence, output_folder/"sequence.fasta", seq_name=self.get_name())
        fm.write_list_to_csv(self.mrf_pos_to_seq_pos, output_folder/"mrf_pos_to_seq_pos.csv")
        fm.write_list_to_csv(self.mrf_pos_to_aln_pos, output_folder/"mrf_pos_to_aln_pos.csv")
        fm.write_list_to_csv(self.aln_pos_to_seq_pos, output_folder/"aln_pos_to_seq_pos.csv")


    def get_insertion_penalties(self):
        return self.insertion_penalties




def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()


    # files
    file_args = parser.add_argument_group('file_args')
    file_args.add_argument('-pf', '--potts_folder', help="Output feature folder", type=pathlib.Path, default=None)
    file_args.add_argument('-aln', '--aln_file', help="Alignment file", type=pathlib.Path)
    file_args.add_argument('-alni', '--aln_with_insertions', help="Alignment file with insertions as lower letters", type=pathlib.Path)
    file_args.add_argument('-ualn', '--unaligned_fasta', help="Unaligned sequences in fasta format", type=pathlib.Path)
    file_args.add_argument('-s', '--sequence_file', help="Sequence file", type=pathlib.Path)
    file_args.add_argument('-p', '--potts_model_file', help="Potts model file", type=pathlib.Path)
    file_args.add_argument('-hhr', '--hhr_file', help="HHblits .hhr output file (needed to find the E-value cutoff)", type=pathlib.Path)

    # hhblits
    hhblits_args = parser.add_argument_group('hhblits_args')
    hhblits_args.add_argument('--hhfilter_threshold', help="HHfilter threshold (default : 80)", type=float, default=80)
    hhblits_args.add_argument('-hhblits', '--call_hhblits', help="Fetch sequences with HHblits", action='store_true', default=False)
    hhblits_args.add_argument('-hd', '--hhblits_database', help="Database path for HHblits call", default=None)
    
    # BLAST
    #parser.add_argument('--nb_sequences_blast', help="Nb sequences fetched by BLAST (default : 100000)", type=int, default=100000)
    #parser.add_argument('--blast_evalue', help="BLAST E-value parameter (default : 1)", type=float, default=1)
    #parser.add_argument('-bxml', '--blast_xml', help="BLAST XML output file (needed to find the E-value cutoff)", type=pathlib.Path)

    # aln files processing
    aln_processing_args = parser.add_argument_group('aln_processing_args')
    aln_processing_args.add_argument('-nofilter', '--dont_filter_alignment', help="Don't filter alignment using HHfilter (default = do)", action='store_true', default=False)
    aln_processing_args.add_argument('-whole', '--use_whole_alignment', help="Use the whole filtered alignment (default : don't : arbitrarily take the first @max_nb_sequences after HHfilter and before trimal)", action='store_true', default=False)
    aln_processing_args.add_argument('-maxnb', '--max_nb_sequences', help="Max. nb sequences in the alignment (if alignment has more sequences that @max_nb_sequences after filtering and before trimming, all sequences after n° @max_nb_sequences will be deleted from the alignment. Default : 1000)", type=int, default=1000)
    aln_processing_args.add_argument('-evcut', '--use_evalue_cutoff', help="Stop taking sequences in the alignment when we reach the elbow of the E-value curve (default : False)", action='store_true', default=False)
    aln_processing_args.add_argument('-notrim', '--dont_trim_alignment', help="Don't trim alignment using trimal (default = do)", action='store_true', default=False)
    aln_processing_args.add_argument('--trimal_gt', help="trimal -gt parameter (default : 0.8)", type=float, default=0.8)
    aln_processing_args.add_argument('--trimal_cons', help="trimal -cons parameter (default : 0)", type=float, default=0)
    aln_processing_args.add_argument('--keep_tmp_files', help="keep temporary files (filtered alignments etc.) (default : false)", action='store_true', default=False)
    aln_processing_args.add_argument('--light', help="keep only Potts model, original sequence and csv files to map positions (default : false)", action='store_true', default=False)

    # limitations
    limitation_args = parser.add_argument_group('limitation_args')
    limitation_args.add_argument('-minnb', '--min_nb_sequences', help="Min. nb sequences in the alignment (if alignment has less than @min_nb_sequences, an exception will be raised and Potts model won't be inferred. Default : 1)", type=int, default=1)
    limitation_args.add_argument('--max_potts_model_length', help="for RAM considerations, won't try to make Potts Object if the train MSA is longer than this (default: 250)", type=int, default=250)


    # CCMpredPy options
    ccmpred_args = parser.add_argument_group('ccmpred_args')
    ccmpred_args.add_argument('--pc_submat', help="CCMpred : Use substitution matrix single pseudocounts instead of uniform (default : False)", default=False, action='store_true')
    ccmpred_args.add_argument('--pc_single_count', help="CCMpred : Specify number of single pseudocounts (default : 1)", default=1)
    ccmpred_args.add_argument('--pc_pair_count', help="CCMpred : Specify number of pair pseudocounts (default : 1)", default=1)
    ccmpred_args.add_argument('--ofn_pll', help="CCMpred : Pseudo-likelihood inference (default : True)", default=True)
    ccmpred_args.add_argument('--ofn_cd', help="CCMpred : Contrastive Divergence inference (default : False)", default=False)
    ccmpred_args.add_argument('--reg_lambda_pair_factor', help="CCMpred : Regularization parameter for pair potentials (L2 regularization with lambda_pair = lambda_pair_factor * scaling) [CCMpred default: 0.2]", default=0.2)
    ccmpred_args.add_argument('--reg_lambda_single', help="CCMpred : Regularization coefficient for single potentials (L2 regularization) [CCMpred default: 10]", default=10)
    ccmpred_args.add_argument('--wt_cutoff', help="CCMpred : Sequence identity threshold. [CCMpred default: 0.8]", default=0.8)


    # Potts model
    potts_model_args = parser.add_argument_group('potts_model_args')
    potts_model_args.add_argument('-noinfer', '--dont_infer_potts_model', help="Don't infer a Potts model (default = do)", action='store_true', default=False)
    potts_model_args.add_argument('--inference_type', help="Inference type (standard : Potts model inferred from an alignment, one_submat : Potts model inferred from a sequence using submatrix pseudocounts, one_hot : one-hot encoding of a sequence -> Potts model) (default : standard)", choices=["standard","one_hot", "one_submat"], default="standard")
    potts_model_args.add_argument('--use_insertion_penalties', help="use insertion penalties (WARNING: EXPERIMENTAL)", action='store_true', default=False)
    potts_model_args.add_argument('--pc_insertions_tau', help="pseudo-count rate for position-specific insertion penalties (WARNING: EXPERIMENTAL)", type=float, default=0)
    
    # Potts model transformation post inference
    post_inference_args = parser.add_argument_group('post_inference_args')
    post_inference_args = parser.add_argument_group('post_inference_args')
    post_inference_args.add_argument('--v_rescaling_function', help="Rescaling function for the v parameters of the Potts model. (default : no rescaling (identity))", default="identity")
    post_inference_args.add_argument('--w_rescaling_function', help="Rescaling function for the w parameters of the Potts model. (default : no rescaling (identity))", default="identity")
    post_inference_args.add_argument('--v_rescaling_tau', help="Tau parameter for rescaling function simulate_uniform_pc_on_v", type=float, default=0.5)
    post_inference_args.add_argument('-nw', '--dont_use_w', help="Speed up computations if we are not interested in w parameters (not recommended)", action='store_true', default=False)
    post_inference_args.add_argument('--insert_null_at_trimmed', help="Insert background parameters at positions trimmed by trimal", action='store_true', default=False)
    post_inference_args.add_argument('--insert_v_star_at_trimmed', help="Insert v parameters at positions trimmed by trimal (computed from frequencies)", action='store_true', default=False)




    args = vars(parser.parse_args(args))

    file_arg_names = ["sequence_file", "aln_file", "aln_with_insertions", "potts_folder", "unaligned_fasta", "potts_model_file", "hhr_file"]
    file_args = {key:args[key] for key in file_arg_names}
    other_args = {key:args[key] for key in args if not key in file_arg_names}

    other_args["filter_alignment"] = not other_args["dont_filter_alignment"]
    other_args["use_less_sequences"] = not other_args["use_whole_alignment"]
    other_args["trim_alignment"] = not other_args["dont_trim_alignment"]
    other_args["infer_potts_model"] = not other_args["dont_infer_potts_model"]
    other_args["use_w"] = not other_args["dont_use_w"]
 
    potts_folder = file_args["potts_folder"]
    if potts_folder is None:
        folder_name = str(uuid.uuid4())
        potts_folder = pathlib.Path(folder_name)
        print("No folder name specified, Potts folder will be created at "+str(potts_folder)) 
    if not potts_folder.is_dir():
        potts_folder.mkdir()

    if other_args["call_hhblits"]:
        potts_object = Potts_Object.from_sequence_with_hhblits(potts_folder, file_args["sequence_file"], **other_args)
    elif file_args["aln_with_insertions"] is not None:
        potts_object = Potts_Object.from_hhblits_files(potts_folder, file_args["aln_with_insertions"], sequence_file=file_args["sequence_file"], potts_model_file=file_args["potts_model_file"], hhr_file=file_args["hhr_file"], **other_args)
    elif file_args["aln_file"] is not None:
        potts_object = Potts_Object.from_aln_file(potts_folder, file_args["aln_file"], sequence_file=file_args["sequence_file"], aln_with_insertions=file_args["aln_with_insertions"], potts_model_file=file_args["potts_model_file"], **other_args)
    elif (file_args["sequence_file"] is not None):
        potts_object = Potts_Object.from_sequence_alone(potts_folder, file_args["sequence_file"], **other_args)
    else:
        raise Exception("TODO")
    
    fm.write_readme(potts_object.folder, **args)

    return potts_object
