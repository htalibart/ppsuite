import unittest
import shutil, tempfile
import pathlib
import pkg_resources

from tests.resources_manager import *

from comutils import files_management as fm
from makepotts.handle_insertions import *
from makepotts.potts_object import *

class Test_Infer_Insertion_Penalties(unittest.TestCase):

    def setUp(self):
        self.potts_folder = pathlib.Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.potts_folder)
        #pass

    def test_count_insertions(self):
       a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"count_insertions_example.a3m")
       ins = count_insertions(a3m_file)
       assert(np.array_equiv(ins,np.array([[0,0,0,0],[0,1,0,0],[0,3,0,0]])))

    def test_infer_insertion_penalties(self):
        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"count_insertions_example.a3m")
        insertion_penalties = infer_insertion_penalties(a3m_file)

    def test_lower_trimmed_columns_insertions(self):
        a3m_file = INSERTION_RESOURCES_FOLDER/'Ac-D.a3m'
        output_file = '/tmp/test_a3m_trim.a3m'
        lower_case_trimmed_columns(a3m_file, output_file, [0,2])
        records = list(SeqIO.parse(str(output_file), 'fasta'))
        assert(str(records[2].seq)=='AcfD')
        os.remove(str(output_file))

    def test_lower_trimmed_columns_insertions_full_object(self):
        a3m_file = INSERTION_RESOURCES_FOLDER/'Ac-D.a3m'
        po = Potts_Object.from_hhblits_files(self.potts_folder, aln_with_insertions=a3m_file, use_insertion_penalties=True, filter_alignment=False, trim_alignment=True, trimal_gt=0.9)
        assert(po.potts_model.ncol==2)
        assert(po.insertion_penalties['open'][1]<po.insertion_penalties['open'][2])
 
    def test_infer_insertion_penalties_with_pseudocounts(self):
        a3m_file = pathlib.Path(INSERTION_RESOURCES_FOLDER/"count_insertions_example.a3m")
        insertion_penalties = infer_insertion_penalties(a3m_file, pc_insertions_tau=0)
        insertion_penalties_with_pseudocounts = infer_insertion_penalties(a3m_file, pc_insertions_tau=0.5)
        print(insertion_penalties)
        print(insertion_penalties_with_pseudocounts)
        assert(insertion_penalties["open"][1]<insertion_penalties_with_pseudocounts["open"][1])
        assert(insertion_penalties["extend"][1]<insertion_penalties_with_pseudocounts["extend"][1])



if __name__=='__main__':
    unittest.main()
