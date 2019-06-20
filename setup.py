from setuptools import setup, find_packages

print(find_packages())

setup(
    name='VAPotts',
    version='0.1dev',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'biopython', 'msgpack'],
    package_data={
        'compotts_solver':['compotts_solver.so'],
        'soeding_reformat':['reformat.pl']
        }
)
