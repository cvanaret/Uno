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



def _glibcxx_version(path):
    """Highest GLIBCXX_3.4.N version provided by a libstdc++ file (N), or None."""
    import re
    try:
        with open(path, 'rb') as file:
            versions = re.findall(rb'GLIBCXX_3\.4\.(\d+)\x00', file.read())
        return max(int(version) for version in versions) if versions else None
    except OSError:
        return None


def _loaded_libstdcxx():
    """Path of the libstdc++ mapped in this process, or None."""
    try:
        with open('/proc/self/maps') as maps:
            for line in maps:
                if '/libstdc++.so' in line:
                    return line.split(None, 5)[-1].strip()
    except OSError:
        pass
    return None


def _preload_bundled_libstdcxx():
    """The Julia runtime of the bundled libkrylov needs a recent libstdc++ (GLIBCXX_3.4.30, GCC 12), and a process can
    hold a single libstdc++.so.6: if an older system libstdc++ is loaded first (e.g. by the unopy extension on
    RHEL 8/9), libkrylov cannot be loaded. Load the bundled (backward-compatible) libstdc++ first, unless the
    libstdc++ that would otherwise be loaded is at least as recent."""
    bundled = os.path.join(basedir, 'libkrylov', 'lib', 'julia', 'libstdc++.so.6')
    if not os.path.isfile(bundled) or _loaded_libstdcxx() is not None:
        return
    bundled_version = _glibcxx_version(bundled)
    if bundled_version is None:
        return
    # the libstdc++ the dynamic linker would load in this environment (same interpreter, LD_LIBRARY_PATH, RPATH)
    import subprocess
    probe = ("import ctypes\n"
             "ctypes.CDLL('libstdc++.so.6')\n"
             "print(next((l.split(None, 5)[-1].strip() for l in open('/proc/self/maps') if '/libstdc++.so' in l), ''))")
    try:
        system = subprocess.run([sys.executable, '-S', '-c', probe], capture_output=True, text=True, timeout=30,
                                check=False).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        system = ''
    system_version = _glibcxx_version(system) if system else None
    if system_version is not None and system_version < bundled_version:
        import ctypes
        try:
            ctypes.CDLL(bundled, mode=ctypes.RTLD_GLOBAL)
        except OSError:
            pass


if sys.platform.startswith('linux'):
    _preload_bundled_libstdcxx()

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
