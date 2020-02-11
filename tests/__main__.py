import unittest
import importlib

test_module_names = ['test_'+n for n in ['makepotts', 'compute_scores', 'call_compotts', 'manage_positions', 'util', 'files_management', 'rescaling']]
test_modules = ['tests.'+tm for tm in test_module_names]
for tm in test_modules:
    globals()[tm] = importlib.import_module(tm)

if __name__=='__main__':
    for tm in test_modules:
        test_module = globals()[tm]
        suite = unittest.TestLoader().loadTestsFromModule(test_module)
        unittest.TextTestRunner(verbosity=3).run(suite)
