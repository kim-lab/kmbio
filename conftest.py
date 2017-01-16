import os
import os.path as op
import atexit

CWD = os.getcwd()
atexit.register(os.chdir, CWD)

os.chdir(op.join(op.dirname(op.abspath(__file__)), 'Tests'))
