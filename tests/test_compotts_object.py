import unittest
import shutil, tempfile
import pathlib
from compotts.compotts_object import *

import pkg_resources
EXAMPLES_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/test_compotts_object/'))

class Test_ComPotts_Object(unittest.TestCase):

    def setUp(self):
        FOLDER = EXAMPLES_FOLDER
        PROTEIN_NAME = "1cc8"
        self.a3m_file = FOLDER/(PROTEIN_NAME+".a3m")
        self.seq_file = FOLDER/(PROTEIN_NAME+".fasta")
        self.output_folder = pathlib.Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.output_folder)

    def test_import_hhblits(self):
        obj = ComPotts_Object(sequence_file=self.seq_file, a3m_file=self.a3m_file, input_folder=self.output_folder, nb_sequences=100)

    def test_to_one_hot(self):
        obj = ComPotts_Object(sequence_file=self.seq_file, input_folder=self.output_folder)

    def test_rescale_compotts_object(self):
        obj = ComPotts_Object(sequence_file=self.seq_file, a3m_file=self.a3m_file, input_folder=self.output_folder, nb_sequences=100, rescaling_function="original_rescaling")


if __name__=='__main__':
    unittest.main()
