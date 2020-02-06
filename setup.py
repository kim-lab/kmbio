import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, setup


def read_md(file):
    with open(file) as fin:
        return fin.read()


PACKAGES = [
    "kmbio",
    "kmbio.PDB",
    "kmbio.PDB.core",
    "kmbio.PDB.io",
    "kmbio.PDB.parsers",
    "kmbio.PDB.tools",
    "kmbio.PDB.tools.QCPSuperimposer",
    "kmbio.SVDSuperimposer",
]

EXTENSIONS = [
    Extension(
        "kmbio.PDB.tools.QCPSuperimposer.qcprotmodule",
        ["kmbio/PDB/tools/QCPSuperimposer/qcprotmodule.c"],
        include_dirs=[np.get_include()],
    ),
    Extension("*", ["kmbio/PDB/parsers/*.pyx"]),
]

setup(
    name="kmbio",
    version="2.0.14",
    author="The Biopython Contributors + KimLab",
    author_email="alexey.strokach@kimlab.org",
    url="https://github.com/kimlaborg/kmbio",
    description="Freely available tools for computational molecular biology.",
    long_description=read_md("README.md"),
    download_url="https://github.com/kimlaborg/kmbortio/release/",
    packages=PACKAGES,
    ext_modules=cythonize(EXTENSIONS),
    package_data={},
)
