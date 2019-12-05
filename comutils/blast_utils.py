import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import pathlib
import subprocess

def call_blast(fastafilename, dbpath, n=100000, blastresultsfile=None, evalue=100):
    if blastresultsfile is None:
        blastresultsfilename = '.'.join(str(fastafilename).split('.')[:-1])+"_blast.xml"
        blastresultsfile = pathlib.Path(blastresultsfilename)

    cline = NcbiblastpCommandline(query=str(fastafilename), db=dbpath, remote=False, out=str(blastresultsfile), outfmt="5", evalue=evalue, max_target_seqs=n)
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
        for blast_record in blast_records:
            if blast_record.descriptions:
                for desc in blast_record.descriptions:
                   results.append(desc.accession)
    return results


def get_blast_record(uid, dbpath):
    tmp_fasta_seq = "temp_seq.fasta"
    cmd = "blastdbcmd -entry "+uid+" -db "+dbpath+" -target_only > "+tmp_fasta_seq
    subprocess.Popen(cmd, shell=True).wait()
    record = list(SeqIO.parse(tmp_fasta_seq, "fasta"))[0]
    os.remove(tmp_fasta_seq)
    return record



def get_full_records_from_blast_xml(blast_xml, dbpath):
    results = get_hits_from_blast_xml(blast_xml)
    hit_records = []
    for hit in results:
        rec = get_blast_record(hit, dbpath)
        hit_records.append(rec)
    return hit_records


def get_fasta_with_full_sequences_from_blast_xml(blast_xml, dbpath, blast_fasta=None):
    if blast_fasta is None:
        blast_fasta_name = '.'.join(str(blast_xml).split('.')[:-1])+".fasta"
        blast_fasta = pathlib.Path(blast_fasta_name)

    hit_records = get_full_records_from_blast_xml(blast_xml, dbpath)
    SeqIO.write(hit_records, blast_fasta, "fasta")
    
    return blast_fasta


def get_blast_xml_and_fasta_output_from_sequence_file(sequence_file, dbpath, blast_fasta=None, blast_xml=None, n=100000, evalue=100):
    blast_xml = call_blast(str(sequence_file), dbpath, n=n, blastresultsfile=blast_xml, evalue=evalue)
    blast_fasta = get_fasta_with_full_sequences_from_blast_xml(blast_xml, dbpath, blast_fasta=blast_fasta)
    return pathlib.Path(blast_xml), pathlib.Path(blast_fasta)
