# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
# flake8: noqa
"""Collection of modules for dealing with biological data in Python.

The Biopython Project is an international association of developers
of freely available Python tools for computational molecular biology.

http://biopython.org
"""
<<<<<<< 7b6395e4f1b4c0fda54a52b277ad68c9d865672b

__version__ = "1.69.4"
=======
__version__ = "1.69.4.dev0"
>>>>>>> Automatically import all modules and move exceptions to a new file

from .exc import *

__all__ = [
    'KDTree',
    'PDB',
    'SVDSuperimposer',
]
from . import *
