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

class Potts_Object:

    def __init__(self):
        pass


    @classmethod
    def from_folder(cls, potts_folder, v_rescaling_function="identity", w_rescaling_function="identity", use_w=True, **kwargs):
        """ instantiate Potts object from Potts folder """
        feature = cls()

        feature.folder = potts_folder
        if not feature.folder.is_dir():
            raise Exception("Feature folder does not exist")

        if (potts_folder/"aln_train.fasta").is_file():
            feature.aln_train = potts_folder/"aln_train.fasta"
        else:
            feature.aln_train = None

        if (potts_folder/"aln_original.fasta").is_file():
            feature.aln_original = potts_folder/"aln_original.fasta"
        else:
            feature.aln_original = None
        

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
        
        return feature




    @classmethod
    def from_files(cls, potts_folder=None, sequence_file=None, potts_model_file=None, aln_file=None, unaligned_fasta=None, fetch_sequences=False, sequences_fetcher='hhblits', database=None, use_evalue_cutoff=False, hhr_file=None, blast_xml=None, filter_alignment=True, hhfilter_threshold=80, use_less_sequences=True, max_nb_sequences=1000, min_nb_sequences=1, trim_alignment=True, trimal_gt=0.8, trimal_cons=0, infer_potts_model=True, inference_type="standard", pc_single_count=1, reg_lambda_pair_factor=0.2, v_rescaling_function="identity", w_rescaling_function="identity", use_w=True, nb_sequences_blast=100000, blast_evalue=1, keep_tmp_files=False, max_potts_model_length=250, insert_null_at_trimmed=False, v_null_is_v0=True, insert_v_star_at_trimmed=False, **kwargs):

        potts_model=None
        if potts_model_file is not None:
            potts_model = Potts_Model.from_msgpack(potts_model_file)

        # ALIGNMENT FOLDER
        if potts_folder is None:
            folder_name = str(uuid.uuid4())
            potts_folder = pathlib.Path(folder_name)
            print("No folder name specified, Potts folder will be created at "+str(potts_folder)) 

        if not potts_folder.is_dir():
            potts_folder.mkdir()

        # FETCH SEQUENCES IF ASKED
        if fetch_sequences:
            if database is None:
                raise Exception("You need to specify a database path (-d option)")
            if sequence_file is None:
                raise Exception("Sequence file missing !")
            if sequences_fetcher=='hhblits':
                aln_file = potts_folder/"aln_original.a3m"
                hhr_file = call_hhblits(sequence_file, aln_file, database, **kwargs)
            elif sequences_fetcher=='blast':
                blast_fasta = potts_folder/"blast.fasta"
                blast_xml = potts_folder/"blast.xml"
                blast_xml, unaligned_fasta = get_blast_xml_and_fasta_output_from_sequence_file(sequence_file, database, blast_fasta=blast_fasta, blast_xml=blast_xml, n=nb_sequences_blast, evalue=blast_evalue)
            else:
                raise Exception(str(sequences_fetcher)+" call not implemented yet. Available options are hhblits or blast")

        # IF EVALUE CUTOFF
        if use_evalue_cutoff:
            if hhr_file is not None:
                cutoff_index = find_hhblits_cutoff_index(hhr_file)
            elif blast_xml is not None:
                cutoff_index = find_blast_cutoff_index(blast_xml)
            else:
                raise Exception("Need a BLAST XML file or .hhr file to find cutoff")


        # ORIGINAL ALIGNMENT, AFTER CUTOFF IF ANY
        if aln_file is not None: # if MSA file is provided
            fm.check_if_file_ok(aln_file)
            if fm.get_format(aln_file)=="fasta":
                aln_original = aln_file
            elif fm.get_format(aln_file)=="a3m": # convert to fasta
                reformat = potts_folder/"reformat.fasta"
                call_reformat(aln_file, reformat)
                aln_original = reformat
            else:
                raise Exception("Unknown format : "+str(fm.get_format(aln_file))+" for "+str(aln_file))
            aln_original_before_clean = aln_original
            aln_original = potts_folder/"aln_original.fasta"
            fm.remove_sequences_with_bad_characters_from_fasta_file_and_upper(aln_original_before_clean, aln_original)
            if use_evalue_cutoff:
                cutoff_fasta = potts_folder/("cutoff_"+str(cutoff_index)+".fasta")
                fm.create_fasta_file_with_less_sequences(aln_original, cutoff_fasta, cutoff_index)
                aln_original = cutoff_fasta
            fm.copy(aln_original, potts_folder/"aln_original.fasta")

        elif unaligned_fasta is not None: # if unaligned sequences are provided
            fm.check_if_file_ok(unaligned_fasta)
            if use_evalue_cutoff:
                cutoff_fasta = potts_folder/("cutoff_"+str(cutoff_index)+".fasta")
                fm.create_fasta_file_with_less_sequences(unaligned_fasta, cutoff_fasta, cutoff_index)
                unaligned_fasta = cutoff_fasta
            if sequence_file is not None:
                fm.add_sequence_to_fasta_file_if_missing(unaligned_fasta, sequence_file)
            clean_unaligned_fasta = potts_folder/"clean_unaligned_sequences.fasta"
            fm.remove_sequences_with_bad_characters_from_fasta_file_and_upper(unaligned_fasta, clean_unaligned_fasta)
            tmp_unaligned_fasta = potts_folder/"clean_unaligned_sequences_bis.fasta"
            fm.copy(clean_unaligned_fasta, tmp_unaligned_fasta)
            fm.create_fasta_file_with_less_sequences(tmp_unaligned_fasta, clean_unaligned_fasta, nb_sequences=10000, fileformat="fasta") # limit to 10000 sequences for MAFFT
            tmp_unaligned_fasta.unlink()
            aln_mafft = potts_folder/"aln_mafft.fasta"
            call_mafft(clean_unaligned_fasta, aln_mafft)
            aln_original = potts_folder/"aln_original.fasta"
            fm.remove_positions_with_gaps_in_first_sequence(aln_mafft, aln_original)
            fm.copy(aln_original, potts_folder/"aln_original.fasta")

        elif inference_type=='one_submat' or inference_type=='one_hot': # if Potts model is inferred from sequence, MSA is set to sequence
            aln_original = sequence_file

        else:
            aln_original = None



        # ORIGINAL MSA TO TRAIN MSA
        aln_train = aln_original

        if (aln_train is not None):
            fm.check_if_file_ok(aln_train)

            if (inference_type=="standard"):
                if filter_alignment: # filter alignment (default keep 80% sequence identity)
                    filtered = potts_folder/"filtered.fasta"
                    call_hhfilter(aln_train, filtered, hhfilter_threshold)
                    aln_train = filtered
                use_less_sequences = use_less_sequences and (not use_evalue_cutoff)
                if use_less_sequences: # if alignment depth fixed arbitrarily (no E-value cutoff)
                    less_fasta = potts_folder/("less_"+str(max_nb_sequences)+".fasta")
                    fm.create_fasta_file_with_less_sequences(aln_train, less_fasta, max_nb_sequences)
                    aln_train = less_fasta

                # trimmed alignment
                if trim_alignment:
                    colnumbering_file = potts_folder/"colnumbering.csv"
                    trimmed_aln = potts_folder/("trim_"+str(int(trimal_gt*100))+".fasta")
                    msa_file_before_trim = aln_train
                    call_trimal(aln_train, trimmed_aln, trimal_gt, trimal_cons, colnumbering_file)
                    aln_train = trimmed_aln
                    mrf_pos_to_aln_pos = fm.get_trimal_ncol(colnumbering_file) 
                else:
                    nb_pos = fm.get_nb_columns_in_alignment(aln_train)
                    mrf_pos_to_aln_pos = [pos for pos in range(nb_pos)]

            else:
                nb_pos = fm.get_nb_columns_in_alignment(aln_train)
                mrf_pos_to_aln_pos = [pos for pos in range(nb_pos)]

            fm.copy(aln_train, potts_folder/"aln_train.fasta")

            if fm.get_nb_sequences_in_fasta_file(aln_train)<min_nb_sequences:
                raise Exception("Less than "+str(min_nb_sequences)+" in the training set : "+str(fm.get_nb_sequences_in_fasta_file(aln_train)))



            # POTTS MODEL
            if (potts_model_file is None) and (infer_potts_model):
                potts_model_file = potts_folder/"potts_model.mrf"

                if fm.get_nb_columns_in_alignment(aln_train)>max_potts_model_length:
                    raise Exception("More than "+str(max_potts_model_length)+" columns in the alignment, won't infer the Potts model.")

                if inference_type=="standard":
                    potts_model = Potts_Model.from_training_set(aln_train, potts_model_file, pc_single_count=pc_single_count, reg_lambda_pair_factor=reg_lambda_pair_factor, **kwargs)

                elif inference_type=="one_submat": # build Potts model from sequence using substitution matrix pseudo-counts on single frequencies
                    potts_model = Potts_Model.from_sequence_file_with_submat(aln_train, filename=potts_model_file, **kwargs)
                
                elif inference_type=="one_hot": # build Potts model from sequence with one-hot encoding
                    potts_model = Potts_Model.from_sequence_file_to_one_hot(aln_train, filename=potts_model_file, **kwargs)
                else:
                    raise Exception("Unknown inference type")


        
        if (potts_model_file is not None):
            if "potts_model" not in locals():
                potts_model = Potts_Model.from_msgpack(potts_model_file)
            potts_model.to_msgpack(potts_model_file)

        if (potts_model_file is not None) and (v_rescaling_function!="identity") and (w_rescaling_function!="identity"):
            if "potts_model" not in locals():
                potts_model = Potts_Model.from_msgpack(potts_model_file)
            potts_model = get_rescaled_potts_model(potts_model, v_rescaling_function, w_rescaling_function, use_w=use_w, **kwargs)
            potts_model.to_msgpack(potts_model_file)



        # IF NO ALN, NO MRF_POS_TO_ALN_POS
        if aln_original is None:
            mrf_pos_to_aln_pos = None


        # SEQUENCE FILE
        if sequence_file is not None:
            fm.copy(sequence_file, potts_folder/"sequence.fasta")
            original_first_seq = fm.get_first_sequence_in_fasta_file(aln_original)
            seq = fm.get_first_sequence_in_fasta_file(sequence_file)
            mrf_pos_to_seq_pos = get_mrf_pos_to_seq_pos(original_first_seq, seq, mrf_pos_to_aln_pos)
            aln_pos_to_seq_pos = get_pos_first_seq_to_second_seq(original_first_seq, seq) 
        elif aln_original is not None: # if no sequence file is provided, sequence is the first sequence of the MSA
            seq = fm.get_first_sequence_in_fasta_file(aln_original)
            seq_name = fm.get_first_sequence_name(aln_original)
            fm.create_seq_fasta(seq, potts_folder/"sequence.fasta", seq_name=seq_name)
            mrf_pos_to_seq_pos = mrf_pos_to_aln_pos
            aln_pos_to_seq_pos = [pos for pos in range(len(seq))]
        else:
            mrf_pos_to_seq_pos = None
            aln_pos_to_seq_pos = None


        # RE-INSERT NULL COLUMNS
        if ((insert_null_at_trimmed) or (insert_v_star_at_trimmed)) and (potts_model is not None):
            if (insert_v_star_at_trimmed):
                potts_model.insert_vi_star_gapped_to_complete_mrf_pos(mrf_pos_to_seq_pos, fm.get_nb_columns_in_alignment(aln_original), msa_file_before_trim)
            elif (insert_null_at_trimmed):
                if v_null_is_v0:
                    v_null = np.tile(get_background_v0(v_rescaling_function, **kwargs), (1,1))
                else:
                    v_null = np.zeros((1,21)) 
                potts_model.insert_null_positions_to_complete_mrf_pos(mrf_pos_to_aln_pos, fm.get_nb_columns_in_alignment(aln_original), v_null=v_null)
            mrf_pos_to_seq_pos = aln_pos_to_seq_pos
            mrf_pos_to_aln_pos = [pos for pos in range(fm.get_nb_columns_in_alignment(aln_original))]
            if (potts_model_file is not None):
                potts_model.to_msgpack(potts_model_file)


        # HANDLE CORRESPONDENCES BETWEEN POTTS MODEL POSITIONS AND SEQ/MSA POSITIONS
        if mrf_pos_to_seq_pos is not None:
            fm.write_list_to_csv(mrf_pos_to_seq_pos, potts_folder/"mrf_pos_to_seq_pos.csv")
        if mrf_pos_to_aln_pos is not None:
            fm.write_list_to_csv(mrf_pos_to_aln_pos, potts_folder/"mrf_pos_to_aln_pos.csv")
        if aln_pos_to_seq_pos is not None:
            fm.write_list_to_csv(aln_pos_to_seq_pos, potts_folder/"aln_pos_to_seq_pos.csv")


        # HANDLE FILES
        if potts_model_file is not None:
            fm.copy(potts_model_file, potts_folder/"potts_model.mrf")
        if not keep_tmp_files:
            for name in ["aln_original.a3m", "reformat.fasta", "filtered.fasta", "less_"+str(max_nb_sequences)+".fasta", "trim_80.fasta", "aln_mafft.fasta", "clean_unaligned_sequences.fasta"]:
                if (potts_folder/name).is_file():
                    (potts_folder/name).unlink()

        return cls.from_folder(potts_folder, v_rescaling_function="identity", w_rescaling_function="identity")



    @classmethod
    def from_merge(cls, potts_folder, objects, aligned_positions_dict, use_less_sequences=False, **kwargs):
        """ create Potts object by merging two Potts objects aligned """
        if not potts_folder.is_dir():
            potts_folder.mkdir()
        aln_original = potts_folder/"aln_original.fasta"
        get_original_msas_aligned_from_aligned_positions(aligned_positions_dict, objects, aln_original)
        return cls.from_files(potts_folder=potts_folder, aln_file=aln_original, use_less_sequences=use_less_sequences, **kwargs)



    def get_name(self):
        if self.sequence_file is not None:
            self.name = fm.get_first_sequence_name(self.sequence_file)
        elif self.aln_original is not None:
            self.name = fm.get_first_sequence_name(self.aln_original)
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
            self.mrf_pos_to_aln_pos = [pos for pos in range(fm.get_nb_columns_in_alignment(aln_original))]


    def to_folder(self, output_folder):
        if not output_folder.is_dir():
            output_folder.mkdir()
        self.potts_model.to_msgpack(output_folder/"potts_model.mrf")
        fm.copy(self.aln_original, output_folder/"aln_original.fasta")
        fm.copy(self.aln_train, output_folder/"aln_train.fasta")
        fm.create_seq_fasta(self.sequence, output_folder/"sequence.fasta", seq_name=self.get_name())
        fm.write_list_to_csv(self.mrf_pos_to_seq_pos, output_folder/"mrf_pos_to_seq_pos.csv")
        fm.write_list_to_csv(self.mrf_pos_to_aln_pos, output_folder/"mrf_pos_to_aln_pos.csv")
        fm.write_list_to_csv(self.aln_pos_to_seq_pos, output_folder/"aln_pos_to_seq_pos.csv")




def main(args=sys.argv[1:]):
    parser = argparse.ArgumentParser()

    # files
    parser.add_argument('-f', '--potts_folder', help="Output feature folder", type=pathlib.Path, default=None)
    parser.add_argument('-aln', '--aln_file', help="Alignment file", type=pathlib.Path)
    parser.add_argument('-ualn', '--unaligned_fasta', help="Unaligned sequences in fasta format", type=pathlib.Path)
    parser.add_argument('-s', '--sequence_file', help="Sequence file", type=pathlib.Path)
    parser.add_argument('-p', '--potts_model_file', help="Potts model file", type=pathlib.Path)

    # fetch sequences ?
    parser.add_argument('-fetch', '--fetch_sequences', help="Fetch sequences in database ? (requires a sequence file) (default : False)", action='store_true', default=False)
    parser.add_argument('-fetcher', '--sequences_fetcher', help="Fetch sequences with...? (hhblits or blast) (default : hhblits)", default='hhblits')
    parser.add_argument('-d', '--database', help="Database path for HHblits or BLAST call", default=None)
    parser.add_argument('--nb_sequences_blast', help="Nb sequences fetched by BLAST (default : 100000)", type=int, default=100000)
    parser.add_argument('--blast_evalue', help="BLAST E-value parameter (default : 1)", type=float, default=1)
    # E-value cutoff
    parser.add_argument('-evcut', '--use_evalue_cutoff', help="Stop taking sequences in the alignment when we reach the elbow of the E-value curve (default : False)", action='store_true', default=False)
    parser.add_argument('-hhr', '--hhr_file', help="HHblits .hhr output file (needed to find the E-value cutoff)", type=pathlib.Path)
    parser.add_argument('-bxml', '--blast_xml', help="BLAST XML output file (needed to find the E-value cutoff)", type=pathlib.Path)

    # Alignment transformation
    parser.add_argument('-nofilter', '--dont_filter_alignment', help="Don't filter alignment using HHfilter (default = do)", action='store_true', default=False)
    parser.add_argument('--hhfilter_threshold', help="HHfilter threshold (default : 80)", type=float, default=80)
    parser.add_argument('-whole', '--use_whole_alignment', help="Use the whole filtered alignment (default : don't : arbitrarily take the first @max_nb_sequences after HHfilter and before trimal)", action='store_true', default=False)
    parser.add_argument('-maxnb', '--max_nb_sequences', help="Max. nb sequences in the alignment (if alignment has more sequences that @max_nb_sequences after filtering and before trimming, all sequences after n° @max_nb_sequences will be deleted from the alignment. Default : 1000)", type=int, default=1000)
    parser.add_argument('-minnb', '--min_nb_sequences', help="Min. nb sequences in the alignment (if alignment has less than @min_nb_sequences, an exception will be raised and Potts model won't be inferred. Default : 1)", type=int, default=1)
    parser.add_argument('-notrim', '--dont_trim_alignment', help="Don't trim alignment using trimal (default = do)", action='store_true', default=False)
    parser.add_argument('--trimal_gt', help="trimal -gt parameter (default : 0.8)", type=float, default=0.8)
    parser.add_argument('--trimal_cons', help="trimal -cons parameter (default : 0)", type=float, default=0)
    parser.add_argument('--keep_tmp_files', help="keep temporary files (filtered alignments etc.) (default : false)", action='store_true', default=False)
    parser.add_argument('--max_potts_model_length', help="for RAM considerations, won't try to make Potts Object if the train MSA is longer than this (default: 250)", type=int, default=250)

    # Potts model
    parser.add_argument('-noinfer', '--dont_infer_potts_model', help="Don't infer a Potts model (default = do)", action='store_true', default=False)
    parser.add_argument('--inference_type', help="Inference type (standard : Potts model inferred from an alignment, one_submat : Potts model inferred from a sequence using submatrix pseudocounts, one_hot : one-hot encoding of a sequence -> Potts model) (default : standard)", default="standard")
    parser.add_argument('--v_rescaling_function', help="Rescaling function for the v parameters of the Potts model. (default : no rescaling (identity))", default="identity")
    parser.add_argument('--w_rescaling_function', help="Rescaling function for the w parameters of the Potts model. (default : no rescaling (identity))", default="identity")
    parser.add_argument('--v_rescaling_tau', help="Tau parameter for rescaling function simulate_uniform_pc_on_v", type=float, default=0.5)
    parser.add_argument('-nw', '--dont_use_w', help="Speed up computations if we are not interested in w parameters (not recommended)", action='store_true', default=False)
    parser.add_argument('--insert_null_at_trimmed', help="Insert background parameters at positions trimmed by trimal", action='store_true', default=False)
    parser.add_argument('--insert_v_star_at_trimmed', help="Insert v parameters at positions trimmed by trimal (computed from frequencies)", action='store_true', default=False)

    # CCMpredPy options
    parser.add_argument('--pc_submat', help="CCMpred : Use substitution matrix single pseudocounts instead of uniform (default : False)", default=False, action='store_true')
    parser.add_argument('--pc_single_count', help="CCMpred : Specify number of single pseudocounts (default : 1)", default=1)
    parser.add_argument('--pc_pair_count', help="CCMpred : Specify number of pair pseudocounts (default : 1)", default=1)
    parser.add_argument('--ofn_pll', help="CCMpred : Pseudo-likelihood inference (default : True)", default=True)
    parser.add_argument('--ofn_cd', help="CCMpred : Contrastive Divergence inference (default : False)", default=False)
    parser.add_argument('--reg_lambda_pair_factor', help="CCMpred : Regularization parameter for pair potentials (L2 regularization with lambda_pair = lambda_pair_factor * scaling) [CCMpred default: 0.2]", default=0.2)



    args = vars(parser.parse_args(args))
    args["filter_alignment"] = not args["dont_filter_alignment"]
    args["use_less_sequences"] = not args["use_whole_alignment"]
    args["trim_alignment"] = not args["dont_trim_alignment"]
    args["infer_potts_model"] = not args["dont_infer_potts_model"]
    args["use_w"] = not args["dont_use_w"]
    cf = Potts_Object.from_files(**args)
    fm.write_readme(cf.folder, **args)
