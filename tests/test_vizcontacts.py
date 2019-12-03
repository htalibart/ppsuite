import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comfeature.comfeature import *
from vizcontacts.contacts_management import *
from vizcontacts.vizpymol import *
from vizcontacts.vizcircos import *

class Test_VizContacts(unittest.TestCase):

    def setUp(self):
        self.feature_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.feature_folder)
        self.comfeature = ComFeature.from_folder(self.feature_folder)

    def tearDown(self):
        #shutil.rmtree(self.feature_folder)
        pass

    def test_get_contact_scores(self):
        contact_scores = get_contact_scores_for_sequence(self.comfeature)

    def test_is_true_contact(self):
        pdb_chain = fm.get_pdb_chain("1cc8", PDB_1CC8)
        assert(is_true_contact((9,67), pdb_chain))

    def test_get_sequence_from_pdb_chain(self):
        pdb_chain = fm.get_pdb_chain("1cc8", PDB_1CC8)
        pdb_sequence = fm.get_sequence_from_pdb_chain(pdb_chain)
        assert(pdb_sequence=="AEIKHYQFNVVMTCSGCSGAVNKVLTKLEPDVSKIDISLEKQLVDVYTTLPYDFILEKIKKTGKEVRSGKQL") 

#    def test_show_n_couplings_pymol(self):
#        contact_scores = get_contact_scores_for_sequence(self.comfeature)
#        pdb_chain = fm.get_pdb_chain("1cc8", PDB_1CC8)
#        real_sequence = fm.get_first_sequence_in_fasta_file(SEQ_1CC8)
#        pdb_contact_scores = translate_dict_to_pdb_pos(contact_scores, pdb_chain, real_sequence)
#        launch_pymol('1cc8', PDB_1CC8)
#        show_n_couplings(50, pdb_contact_scores, PDB_1CC8, '1cc8')

#    def test_show_pymol(self):
#        shutil.copy(PDB_1CC8, self.feature_folder)
#        show_predicted_contacts_with_pymol(self.feature_folder, "1cc8", chain_id='A', coupling_sep_min=3)

#    def test_create_circos(self):
#        circos_output_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
#        coupling_dicts_for_sequence_indexed_by_colors = {"red":{(1,50):0.6}, "blue":{(5,7):0.9}}
#        sequence = self.comfeature.sequence
#        create_circos(circos_output_folder, coupling_dicts_for_sequence_indexed_by_colors, sequence)


    def test_create_circos_from_comfeature_and_pdb_chain(self):
        pdb_chain = fm.get_pdb_chain("1cc8", PDB_1CC8)
        create_circos_from_comfeature_and_pdb_chain(self.comfeature, pdb_chain)

if __name__=='__main__':
    unittest.main()
