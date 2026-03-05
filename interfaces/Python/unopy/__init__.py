import os
import sys

# load bundled DLLs from unopy.libs
basedir = os.path.dirname(__file__)
subdir = os.path.join(basedir, '..', 'unopy.libs')
if os.name == 'nt':
    os.add_dll_directory(os.path.abspath(subdir))
elif sys.platform == 'cygwin':
    os.environ['PATH'] = os.pathsep.join((os.environ['PATH'], basedir, subdir))

from .unopy import *  # re-export the .pyd contents