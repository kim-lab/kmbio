"""py.test configuration file

This file is required because tests in the Tests folder expect
the PDB folder to be relative to the test file.
"""
import atexit
import os
import os.path as op

CWD = os.getcwd()
TEST_DIR = op.dirname(op.abspath(__file__))
if CWD != TEST_DIR:
    os.chdir(TEST_DIR)
    atexit.register(os.chdir, CWD)
