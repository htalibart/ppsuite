from setuptools import setup, find_packages

print(find_packages())

setup(
    name='VAPotts',
    version='0.1dev',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'biopython', 'msgpack==0.6.1', 'msgpack-python==0.5.6', 'scipy', 'matplotlib', 'seaborn', 'sklearn', 'kneebow'],
    package_data={
        'compotts':['compotts_solver.so'],
        'makepotts':['P2Pmat14Weighted_prob.csv', 'P2Pmat14WeightedInterIntra=_prob_cond.csv', 'P2Pmat14WeightedInter=_prob_cond.csv', 'P2Pmat14Inter=_prob_cond.csv'],
        },
    entry_points={
        'console_scripts':['compotts = compotts.__main__:main', 'vizpotts = vizpotts.__main__:main', 'makemsa = makemsa.__main__:main', 'makepotts = makepotts.potts_object:main', 'vizpymol = vizcontacts.vizpymol:main', 'vizcircos = vizcontacts.vizcircos:main', 'vizmap = vizcontacts.vizmap:main']
        },
    test_suite="tests",
)
