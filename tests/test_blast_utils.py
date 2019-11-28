import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from basic_modules.blast_utils import *


class Test_BLAST_Utils(unittest.TestCase):

    def test_get_hits_from_blast_xml(self):
        records = get_hits_from_blast_xml(BLAST_XML)
        assert(len(records)>0)

    def test_get_full_records_from_blast_xml(self):
        full_records = get_full_records_from_blast_xml(BLAST_XML)
        assert(len(full_records)>0)

    def test_get_fasta_with_full_sequences_from_blast_xml(self):
        blast_fasta = get_fasta_with_full_sequences_from_blast_xml(BLAST_XML)
        assert(blast_fasta.exists())

if __name__=='__main__':
    unittest.main()
