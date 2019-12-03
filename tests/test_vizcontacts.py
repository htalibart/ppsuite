import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comfeature.comfeature import *
from vizcontacts.contacts_management import *
from vizcontacts.vizpymol import *

class Test_VizContacts(unittest.TestCase):

    def setUp(self):
        self.feature_folder = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        shutil.copytree(FEATURE_FOLDER, self.feature_folder)
        self.comfeature = ComFeature.from_folder(self.feature_folder)

    def tearDown(self):
        shutil.rmtree(self.feature_folder)

    def test_get_contact_scores(self):
        contact_scores = get_contact_scores_for_sequence(self.comfeature)

    def test_is_true_contact(self):
        pdb_chain = fm.get_pdb_chain("1cc8", PDB_1CC8)
        assert(is_true_contact((9,67), pdb_chain))

    def test_show_n_couplings_pymol(self):
        contact_scores = get_contact_scores_for_sequence(self.comfeature)
        launch_pymol('1cc8', PDB_1CC8)
        show_n_couplings(50, contact_scores, PDB_1CC8, '1cc8')

if __name__=='__main__':
    unittest.main()
