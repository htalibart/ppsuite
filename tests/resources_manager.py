import pathlib
import pkg_resources

RESOURCES_1CC8_EVERYTHING_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/1cc8_everything/'))
SEQ_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.fasta"
A3M_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.a3m"
ALN_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_trim_80.fasta"
SMALL_ALN_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_small_aln.fasta"
LARGER_SMALL_ALN_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_larger_small_aln.fasta"
LARGER_SMALL_ALN_1CC8_TRIM_80 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_larger_small_aln_trim_80.fasta"
MRF_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_standard.mrf"
PDB_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.pdb"

FEATURE_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/feature_folder/'))

FAKE_SEQS_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/fake_sequences/'))

BLAST_XML = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/blast.xml'))
BLAST_FASTA = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/blast.fasta'))

FOLDER_5JZR = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/5jzr/'))
A3M_5JZR = FOLDER_5JZR/"Q9AZ42.a3m"
HHR_5JZR = FOLDER_5JZR/"Q9AZ42.hhr"
SEQ_5JZR = FOLDER_5JZR/"Q9AZ42.fasta"
PDB_5JZR = FOLDER_5JZR/"5jzr.pdb"

TEST_TS_MSA = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/ts_msa.fasta'))
TEST_TS_SEQ = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/ts_seq.fasta'))
TEST_TS_SEQ_TRIM = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/ts_seq_trim.fasta'))
