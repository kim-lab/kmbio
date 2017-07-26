import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup

PACKAGES = [
    'kmbio',
    'kmbio.PDB',
    'kmbio.PDB.core',
    'kmbio.PDB.parsers',
    'kmbio.PDB.tools',
    'kmbio.PDB.tools.QCPSuperimposer',
    'kmbio.SVDSuperimposer',
    'kmbio.KDTree',
]

EXTENSIONS = [
    Extension(
        'kmbio.KDTree._CKDTree',
        ["kmbio/KDTree/KDTree.c", "kmbio/KDTree/KDTreemodule.c"],
        include_dirs=[np.get_include()]),
    Extension(
        'kmbio.PDB.tools.QCPSuperimposer.qcprotmodule',
        ["kmbio/PDB/tools/QCPSuperimposer/qcprotmodule.c"],
        include_dirs=[np.get_include()]),
    Extension(
        '*',
        ["kmbio/PDB/parsers/*.pyx"]),
]


setup(
    name='kmbio',
    version="1.69.6.dev0",
    author='The Biopython Contributors + KimLab',
    author_email='alexey.strokach@kimlab.org',
    url='https://github.com/kimlaborg/kmbio',
    description="Freely available tools for computational molecular biology.",
    download_url='https://github.com/kimlaborg/kmbio/release/',
    packages=PACKAGES,
    ext_modules=cythonize(EXTENSIONS),
    package_data={},
)
