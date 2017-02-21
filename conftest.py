import atexit
import os
import os.path as op
import sys

CWD = os.getcwd()
atexit.register(os.chdir, CWD)

# Use the installed version of the package instead of the current directory
sys.path.remove(CWD)

os.chdir(op.join(op.dirname(op.abspath(__file__)), 'Tests'))
