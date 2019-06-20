from setuptools import setup, find_packages

setup(
    name='VAPotts',
    version='0.1dev',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'biopython'], # TODO
    long_description=open('README.md').read(),
)
