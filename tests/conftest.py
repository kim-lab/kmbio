"""py.test configuration file

This file is required because tests in the Tests folder expect
the PDB folder to be relative to the test file.
"""
import atexit
import logging
import os
import os.path as op

import pytest

logger = logging.getLogger(__name__)

CWD = os.getcwd()
TEST_DIR = op.dirname(op.abspath(__file__))
if CWD != TEST_DIR:
    os.chdir(TEST_DIR)
    atexit.register(os.chdir, CWD)


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
