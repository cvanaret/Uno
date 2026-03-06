import os
import sys

# load bundled DLLs from unopy.libs
basedir = os.path.dirname(__file__)
subdir = os.path.join(basedir, '..', 'unopy.libs')
if os.name == 'nt':
    # Prepend to PATH so our DLLs are found before any system MinGW ones
    os.environ['PATH'] = subdir + os.pathsep + os.environ.get('PATH', '')
    os.add_dll_directory(subdir)
elif sys.platform == 'cygwin':
    os.environ['PATH'] = os.pathsep.join((os.environ['PATH'], basedir, subdir))

from .unopy import *