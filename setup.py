from setuptools import setup, find_packages

print(find_packages())

setup(
    name='VAPotts',
    version='0.1dev',
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'biopython', 'msgpack'],
    package_data={
        'compotts':['compotts_solver.so'],
        },
    entry_points={
        'console_scripts':['compotts = compotts.__main__:main']
        },
    test_suite="tests",
)
