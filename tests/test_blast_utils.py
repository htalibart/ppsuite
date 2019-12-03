# TODO CHANGER LES FICHIERS DE TEST
import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils.blast_utils import *


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

    def test_call_blast(self):
        blast_xml = pathlib.Path('/tmp/'+next(tempfile._get_candidate_names()))
        nb_seq = 10
        call_blast(SEQ_1CC8, "nr", n=nb_seq, blastresultsfile=blast_xml, evalue=0.1, remote=True)
        records = get_hits_from_blast_xml(blast_xml)
        print(len(records))
        assert(len(records)==nb_seq)

if __name__=='__main__':
    unittest.main()
