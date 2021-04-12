from setuptools import setup, find_packages

print(find_packages())

setup(
    name='PPSuite',
    version='1.0',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'biopython', 'msgpack==0.6.1', 'msgpack-python==0.5.6', 'scipy', 'matplotlib', 'seaborn', 'sklearn', 'kneebow'],
    package_data={
        'ppalign':['ppalign_solver.so'],
        'tests':['examples/*', 'examples/*/*'],
        'makepotts':['infer_insertion_penalties.jl'],
        },
    entry_points={
        'console_scripts':['ppalign = ppalign.__main__:main', 'vizpotts = vizpotts.__main__:main', 'makepotts = makepotts.potts_object:main', 'vizpymol = vizcontacts.vizpymol:main', 'vizcircos = vizcontacts.vizcircos:main', 'vizmap = vizcontacts.vizmap:main', 'fasta2csv = comutils.fasta2indices:main']
        },
    test_suite="tests",
)
