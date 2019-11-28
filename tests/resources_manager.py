import pathlib
import pkg_resources

RESOURCES_1CC8_EVERYTHING_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/1cc8_everything/'))
SEQ_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.fasta"
A3M_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8.a3m"
MRF_1CC8 = RESOURCES_1CC8_EVERYTHING_FOLDER/"1cc8_standard.mrf"
FAKE_SEQS_FOLDER = pathlib.Path(pkg_resources.resource_filename(__name__,'examples/fake_sequences/'))
