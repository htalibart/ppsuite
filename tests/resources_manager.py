import pathlib
import pkg_resources

RESOURCES_1CC8_EVERYTHING_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/1cc8_everything/'))
SEQ_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.fasta"
A3M_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.a3m"
MRF_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_standard.mrf"

FAKE_SEQS_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/fake_sequences/'))

BLAST_XML = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/example_blast.xml'))

TEST_TS_MSA = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/ts_msa.fasta'))
TEST_TS_SEQ = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/ts_seq.fasta'))
TEST_TS_SEQ_TRIM = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/ts_seq_trim.fasta'))
