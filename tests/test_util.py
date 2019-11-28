import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from basic_modules.util import *
from basic_modules import files_management as fm


class Test_Util(unittest.TestCase):

    def test_get_trimmed_sequence_for_msa(self):
        seq = fm.get_first_sequence_in_fasta_file(TEST_TS_SEQ)
        seq_trim = get_trimmed_sequence_for_msa(TEST_TS_MSA, seq)
        real_seq_trim = fm.get_first_sequence_in_fasta_file(TEST_TS_SEQ_TRIM)
        assert(seq_trim==real_seq_trim)

if __name__=='__main__':
    unittest.main()
