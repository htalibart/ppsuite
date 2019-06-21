import unittest
import test_compotts_object, test_call_compotts_simple, test_call_compotts_one_hot, test_call_compotts_hhblits_object, test_call_compotts_one_seq, test_manage_positions

if __name__=='__main__':
    test_modules = [test_compotts_object, test_call_compotts_simple, test_call_compotts_one_hot, test_call_compotts_hhblits_object, test_call_compotts_one_seq, test_manage_positions]
    for test_module in test_modules:
        suite = unittest.TestLoader().loadTestsFromModule(test_module)
        unittest.TextTestRunner(verbosity=3).run(suite)
