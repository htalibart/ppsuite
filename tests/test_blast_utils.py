import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils.blast_utils import *
from comutils import files_management as fm

class Test_BLAST_Utils(unittest.TestCase):
    def setUp(self):
        self.db_path = '~/data/blastdb/swissprot'

    def test_get_hits_from_blast_xml(self):
        records = get_hits_from_blast_xml(BLAST_XML)
        assert(len(records)>0)

    def test_get_blast_xml_and_fasta_output_from_sequence_file(self):
        blast_xml = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        blast_fasta = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        nb_sequences=100
        blast_xml, blast_fasta = get_blast_xml_and_fasta_output_from_sequence_file(SEQ_1CC8, self.db_path, blast_fasta=blast_fasta, blast_xml=blast_xml, n=nb_sequences, evalue=100)
        assert(blast_xml.exists())
        assert(blast_fasta.exists())
        assert(fm.get_nb_sequences_in_fasta_file(blast_fasta)==nb_sequences)


if __name__=='__main__':
    unittest.main()
