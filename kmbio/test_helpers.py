import logging
import random

import pytest

logger = logging.getLogger(__name__)

# Constants
PDB_IDS = ["4dkl", "1arr", "1dvf", "3mbp"]

#: There are no precalculated mmCIF biounit structures for PDBs which have no biounits.
#: Conversely, PDBs which make up a large biounit often only have mmCIF structure.
MISSING = [
    # (pdb_id, pdb_type, biounit)
    ("1arr", "pdb", 1),
    ("4p6f", "pdb", 0),
    ("4p6f", "pdb", 1),
]

DIFFICULT = ["4p6f"]

LOCAL_REMOTE_MISMATCH = [("4dkl", "pdb"), ("4dkl", "pdb")]

# PDBs that cause errors
ATOM_DEFINED_TWICE_PDBS = ["2q3u", "2kax", "1wcn", "1wco", "2dii", "2eya"]

NO_RESNAME_ATTRIBUTE_PDBS = ["1q3l", "4d1e", "1cty", "4pru", "1ctz", "2h9p"]

NUM_TESTS_MAX = 4


def set_random_seed():
    random.seed(42)


def random_subset(lst, nmax=NUM_TESTS_MAX):
    random.shuffle(lst)
    return lst[:nmax]


def parametrize(arg_string, arg_list):
    """
    Args:
        arg_string: Comma-separated string of arguments (e.g. 'pdb_id, pdb_type').
        arg_list: List of arguments or argument dictionaries.
    """
    logger.info("arg_string: %s", arg_string)
    logger.info("arg_list: %s", arg_list)
    if "," in arg_string:
        keys = arg_string.replace(" ", "").split(",")
        args = [tuple(r[k] for k in keys) for r in arg_list]
    else:
        key = arg_string.replace(" ", "")
        if isinstance(arg_list[0], dict) and key in arg_list[0]:
            args = [r[key] for r in arg_list]
        else:
            args = arg_list
    logger.info("args: %s", args)
    return pytest.mark.parametrize(arg_string, args)
