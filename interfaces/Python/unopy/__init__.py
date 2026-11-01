import os
import sys

# Bundled DLLs ship in two places:
#   unopy/lib    - libkrylov + the Julia runtime (installed by CMake)
#   unopy.libs   - libraries vendored by delvewheel (e.g. the Fortran runtime)
basedir = os.path.dirname(__file__)
libdir = os.path.join(basedir, 'lib')
subdir = os.path.join(basedir, '..', 'unopy.libs')

if os.name == 'nt':
    # Prepend to PATH so our DLLs win over any system MinGW ones. libdir is
    # prepended last, so it lands first on PATH — Julia's libstdc++ should take
    # precedence over the one delvewheel vendored from msys2.
    for d in (subdir, libdir):
        if os.path.isdir(d):
            os.environ['PATH'] = d + os.pathsep + os.environ.get('PATH', '')
            os.add_dll_directory(d)
elif sys.platform == 'cygwin':
    os.environ['PATH'] = os.pathsep.join((os.environ['PATH'], basedir, libdir, subdir))

from .unopy import *