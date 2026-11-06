import os
import sys

# load bundled DLLs from unopy.libs
basedir = os.path.dirname(__file__)
subdir = os.path.join(basedir, '..', 'unopy.libs')
if os.name == 'nt' and os.path.isdir(subdir):
    # Prepend to PATH so our DLLs are found before any system MinGW ones
    os.environ['PATH'] = subdir + os.pathsep + os.environ.get('PATH', '')
    os.add_dll_directory(subdir)
elif sys.platform == 'cygwin':
    os.environ['PATH'] = os.pathsep.join((os.environ['PATH'], basedir, subdir))

from .unopy import *
from .unopy import _set_default_libkrylov_path


def bundled_libkrylov_path():
    """Path to the libkrylov library bundled with unopy (Python wheels), or None if there is none."""
    if sys.platform == 'win32':
        candidate = os.path.join(basedir, 'libkrylov', 'bin', 'libkrylov.dll')
    elif sys.platform == 'darwin':
        candidate = os.path.join(basedir, 'libkrylov', 'lib', 'libkrylov.dylib')
    else:
        candidate = os.path.join(basedir, 'libkrylov', 'lib', 'libkrylov.so')
    return candidate if os.path.isfile(candidate) else None


# libkrylov is loaded lazily, only if needed. Resolution order: option libkrylov_path > UNO_KRYLOV_LIBRARY environment
# variable > bundled library > platform default
if bundled_libkrylov_path() is not None:
    _set_default_libkrylov_path(bundled_libkrylov_path())
