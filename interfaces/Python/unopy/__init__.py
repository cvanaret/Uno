import os
import sys

# load bundled DLLs from unopy.libs
_libs_dir = os.path.join(os.path.dirname(__file__), '..', 'unopy.libs')
if os.path.isdir(_libs_dir):
    if sys.version_info >= (3, 8):
        os.add_dll_directory(os.path.abspath(_libs_dir))
    else:
        os.environ['PATH'] = os.path.abspath(_libs_dir) + os.pathsep + os.environ.get('PATH', '')

from .unopy import *  # re-export the .pyd contents