import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from makepotts.potts_object import *
from comutils import files_management as fm

class Test_MakePotts(unittest.TestCase):

    def setUp(self):
        self.potts_folder = pathlib.Path(tempfile.mkdtemp())

    def tearDown(self):
        #shutil.rmtree(self.potts_folder)
        pass

    def test_makepotts_with_potts_model_only(self):
       potts_model_file = pathlib.Path(MRF_1CC8)
       cf = Potts_Object.from_potts_model(self.potts_folder, potts_model_file=potts_model_file)

    def test_makepotts_with_aln_only_simple(self):
       aln_file = pathlib.Path(LARGER_SMALL_ALN_1CC8)
       cf = Potts_Object.from_aln_file(self.potts_folder, aln_file=aln_file, trim_alignment=False)
       assert(cf.potts_model.ncol==fm.get_nb_columns_in_alignment(aln_file))

    def test_pos_trim_aln(self):
       aln_file = pathlib.Path(LARGER_SMALL_ALN_1CC8)
       aln_file_trim_80 = pathlib.Path(LARGER_SMALL_ALN_1CC8_TRIM_80)
       cf = Potts_Object.from_aln_file(self.potts_folder, aln_file=aln_file, trim_alignment=True, trimal_gt=0.8, trimal_cons=0)
       assert(cf.potts_model.ncol==fm.get_nb_columns_in_alignment(aln_file_trim_80))
       assert(len(cf.mrf_pos_to_aln_pos)==fm.get_nb_columns_in_alignment(aln_file_trim_80))
       assert(len(cf.mrf_pos_to_seq_pos)==fm.get_nb_columns_in_alignment(aln_file_trim_80))

    def test_correctly_trimmed(self):
       aln_file = pathlib.Path(LARGER_SMALL_ALN_1CC8)
       seq_file = pathlib.Path(SEQ_1CC8)
       cf = Potts_Object.from_aln_file(self.potts_folder, aln_file=aln_file, trim_alignment=True, trimal_gt=0.8, trimal_cons=0, sequence_file=seq_file)
       sequence = fm.get_first_sequence_in_fasta_file(seq_file)
       sequence_mrf = ''.join(sequence[ind] for ind in cf.mrf_pos_to_seq_pos)
       aln_file_trim_80 = pathlib.Path(LARGER_SMALL_ALN_1CC8_TRIM_80)
       sequence_trim = fm.get_first_sequence_in_fasta_file(aln_file_trim_80)
       assert(sequence_mrf==sequence_trim)

    def test_cutoff_hhr(self):
       a3m_file = pathlib.Path(A3M_5JZR)
       hhr_file = pathlib.Path(HHR_5JZR)
       cf = Potts_Object.from_hhblits_files(self.potts_folder, aln_with_insertions=a3m_file, hhr_file=hhr_file, use_evalue_cutoff=True, infer_potts_model=False)
       nb_sequences = fm.get_nb_sequences_in_fasta_file(cf.aln_train)
       assert(nb_sequences==7366)

    def test_less_sequences(self):
       cf = Potts_Object.from_aln_file(self.potts_folder, aln_file=ALN_1CC8, trim_alignment=False, use_less_sequences=True, max_nb_sequences=1000, infer_potts_model=False)
       nb_sequences = fm.get_nb_sequences_in_fasta_file(cf.aln_train)
       assert(nb_sequences==1000)

    def test_nb_min_sequences(self):
       with self.assertRaises(Exception) as context:
           cf = Potts_Object.from_aln_file(self.potts_folder, aln_file=ALN_1CC8, trim_alignment=False, use_less_sequences=True, max_nb_sequences=1000, min_nb_sequences=2000, infer_potts_model=False)

    def test_insert_null_position(self):
       pos=5
       potts_model = Potts_Model.from_msgpack(pathlib.Path(MRF_1CC8))
       L = potts_model.ncol
       potts_model.insert_null_position_at(pos)
       assert(potts_model.ncol==L+1)
       assert(potts_model.v.shape==(L+1,21))
       assert(potts_model.w.shape==(L+1,L+1,21,21))
       assert(potts_model.v[pos][5]==0)
       assert(potts_model.w[0,pos][1,7]==0)

    def test_insert_null_positions(self):
       self.potts_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
       shutil.copytree(FEATURE_FOLDER, self.potts_folder)
       cf = Potts_Object.from_folder(self.potts_folder)
       cf.potts_model.insert_null_positions_to_complete_mrf_pos(cf.mrf_pos_to_seq_pos, len(cf.sequence))

    def test_makepotts_with_null_positions(self):
       cf = Potts_Object.from_hhblits_files(self.potts_folder, A3M_1CC8, insert_null_at_trimmed=True)
       assert(cf.potts_model.ncol==len(cf.sequence))

    #def test_from_file_calling_hhblits_and_evalue_cutoff(self):
    #    cf = Potts_Object.from_sequence_with_hhblits(self.potts_folder, SEQ_1CC8, '~/data/uniclust30_2018_08/uniclust30_2018_08', use_evalue_cutoff=True) 


    def test_cutoff_blast(self):
       cf = Potts_Object.from_blast_files(self.potts_folder, SEQ_1CC8, unaligned_fasta=BLAST_FASTA, blast_xml=BLAST_XML, use_evalue_cutoff=True, infer_potts_model=False, filter_alignment=False)
       nb_sequences = fm.get_nb_sequences_in_fasta_file(cf.aln_train)
       assert(nb_sequences==17) # because clean

    def test_one_submat(self):
       cf = Potts_Object.from_sequence_alone(self.potts_folder, SEQ_1CC8, inference_type="one_submat")
       assert(cf.potts_model_file.is_file())

    def test_one_hot(self):
       cf = Potts_Object.from_sequence_alone(self.potts_folder, SEQ_1CC8, inference_type="one_hot")
       assert(cf.potts_model_file.is_file())


    def test_insertion_penalties(self):
       a3m_file = pathlib.Path(A3M_5JZR)
       po = Potts_Object.from_hhblits_files(self.potts_folder, aln_with_insertions=a3m_file, use_insertion_penalties=True)
       assert(po.insertion_penalties is not None)


    def test_lower_trimmed_columns_insertions(self):
        a3m_file = INSERTION_RESOURCES_FOLDER/'Ac-D.a3m'
        output_file = '/tmp/test_a3m_trim.a3m'
        lower_case_trimmed_columns(a3m_file, output_file, [0,2])
        records = list(SeqIO.parse(str(output_file), 'fasta'))
        assert(str(records[2].seq)=='AcfD')
        os.remove(str(output_file))

    def test_lower_trimmed_columns_insertions_full_object(self):
        a3m_file = INSERTION_RESOURCES_FOLDER/'Ac-D.a3m'
        po = Potts_Object.from_hhblits_files(self.potts_folder, aln_with_insertions=a3m_file, use_insertion_penalties=True, filter_alignment=False, trim_alignment=True, trimal_gt=0.9)
        assert(po.potts_model.ncol==2)
        assert(po.insertion_penalties['open'][1]<po.insertion_penalties['open'][2])
 
        
    def test_zero_fill(self):
        L=3
        p = Potts_Model.zero_fill(L)
        assert(p.ncol==L)
        assert(p.get_v_norm()==0)
        assert(p.get_w_norm()==0)



    def test_infer_with_adabmdca(self):
       aln_file = pathlib.Path(LARGER_SMALL_ALN_1CC8)
       po_ada = Potts_Object.from_aln_file(self.potts_folder, aln_file=aln_file, inference_method='adabmDCA')
       po_ccm = Potts_Object.from_aln_file(self.potts_folder, aln_file=aln_file, inference_method='CCMpredPy')
       assert(po_ada.potts_model.ncol==po_ccm.potts_model.ncol)

if __name__=='__main__':
    unittest.main()
