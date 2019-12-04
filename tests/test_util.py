import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils.util import *
from comutils import files_management as fm


class Test_Util(unittest.TestCase):

    def test_get_trimmed_sequence_for_msa(self):
        seq = fm.get_first_sequence_in_fasta_file(TEST_TS_SEQ)
        seq_trim = get_trimmed_sequence_for_msa(TEST_TS_MSA, seq)
        real_seq_trim = fm.get_first_sequence_in_fasta_file(TEST_TS_SEQ_TRIM)
        assert(seq_trim==real_seq_trim)

    def test_get_pos_dict_first_seq_to_second_seq(self): #d[pos_in_first_seq] = pos_in_second_seq 
        first_seq = "ACDEFGH"
        begin_index = 1
        second_seq = first_seq[begin_index:]
        d = get_pos_first_seq_to_second_seq(first_seq, second_seq)
        for pos_in_first_seq in range(0,begin_index):
            pos_in_second_seq = d[pos_in_first_seq]
            assert(pos_in_second_seq is None)
        for pos_in_first_seq in range(begin_index,len(first_seq)):
            pos_in_second_seq = d[pos_in_first_seq]
            assert(first_seq[pos_in_first_seq]==second_seq[pos_in_second_seq])


    def test_get_pos_dict_first_seq_to_second_seq_1cc8(self): 
        real_seq = fm.get_first_sequence_in_fasta_file(SEQ_1CC8)
        pdb_seq = fm.get_sequence_from_pdb_chain(fm.get_pdb_chain('1cc8', PDB_1CC8, chain_id='A'))
        d = get_pos_first_seq_to_second_seq(real_seq, pdb_seq)
        for pos_in_real_seq in range(len(real_seq)):
            if d[pos_in_real_seq] is not None:
                assert(pdb_seq[d[pos_in_real_seq]]==real_seq[pos_in_real_seq])


    def test_get_pos_dict_first_seq_to_second_seq_5jzr(self): 
        real_seq = fm.get_first_sequence_in_fasta_file(SEQ_5JZR)
        pdb_seq = fm.get_sequence_from_pdb_chain(fm.get_pdb_chain('5jzr', PDB_5JZR, chain_id='A'))
        d = get_pos_first_seq_to_second_seq(real_seq, pdb_seq)
        for pos_in_real_seq in range(len(real_seq)):
            if d[pos_in_real_seq] is not None:
                assert(pdb_seq[d[pos_in_real_seq]]==real_seq[pos_in_real_seq])


if __name__=='__main__':
    unittest.main()
