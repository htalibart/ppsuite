from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import pathlib

from comutils.db_utils import *

def call_blast(fastafilename, dbpath, n=100000, blastresultsfile=None, evalue=1000000000):
    # TODO test
    if blastresultsfile is None:
        blastresultsfilename = '.'.join(str(fastafilename).split('.')[:-1])+"_blast.xml"
        blastresultsfile = pathlib.Path(blastresultsfilename)

    cline = NcbiblastpCommandline(query=fastafilename, db=dbpath, remote=False, out=blastresultsfilename, outfmt="5", evalue=evalue, max_target_seqs=n)
    print('Running:', cline)
    stdout, stderr = cline()
    if stdout:
        print('stdout:')
        print(stdout)
    if stderr:
        print('stderr')
        print(stderr)
    
    return blastresultsfile


def get_hits_from_blast_xml(blast_xml):
    results = []
    with open(blast_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for recordcount, blast_record in enumerate(blast_records):
            if blast_record.descriptions:
                for desc in blast_record.descriptions:
                    results.append(desc.title.split(' ')[0])
    return results


def get_full_records_from_blast_xml(blast_xml):
    sprot = SwissprotAccessor(sprot_filepath=None)
    results = get_hits_from_blast_xml(blast_xml)
    hit_records = []
    for hit in results:
        uid = hit.split('|')[1]
        try:
            rec = sprot.get_record(uid)
            hit_records.append(rec)
        except Exception as e:
            print(hit, e)
    return hit_records


def get_fasta_with_full_sequences_from_blast_xml(blast_xml, blast_fasta=None):
    if blast_fasta is None:
        blast_fasta_name = '.'.join(str(blast_xml).split('.')[:-1])+".fasta"
        blast_fasta = pathlib.Path(blast_fasta_name)

    hit_records = get_full_records_from_blast_xml(blast_xml)
    SeqIO.write(hit_records, blast_fasta, "fasta")
    
    return blast_fasta
