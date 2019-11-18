from setuptools import setup, find_packages

print(find_packages())

setup(
    name='VAPotts',
    version='0.1dev',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'biopython', 'msgpack', 'scipy', 'matplotlib', 'seaborn'],
    package_data={
        'compotts':['compotts_solver.so'],
        },
    entry_points={
        'console_scripts':['compotts = compotts.__main__:main', 'vizpotts = vizpotts.__main__:main', 'filepotts = compotts.align_msas:main', 'multi_compotts = multi_compotts.__main__:main', 'comfeature = compotts.compotts_object:main']
        },
    test_suite="tests",
)
