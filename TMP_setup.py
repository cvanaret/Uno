"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""
import os
import platform
import re
import subprocess
import sys

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

def get_machine():
    # x86_64 / AMD64 / arm64 
    machine = platform.machine().lower()
    if machine == "x86_64" or machine == "amd64":
        return "x86_64"
    elif machine == "arm64":
        return "aarch64"
    else:
        raise RuntimeError("The machine " + machine + " is not known.")

def get_system():
    # linux / win32 / darwin
    system = sys.platform
    if system == "linux":
        return "linux-gnu"
    elif system == "win32":
        return "w64-mingw32"
    elif system == "darwin":
        return "apple-darwin"
    else:
        raise RuntimeError("The system " + system + " is not known.")

# The main interface is through CMake
class CMakeBuild(build_ext):
    def build_extension(self, ext):
        # paths
        current_directory = os.getcwd()
        
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # required for for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep
        
        ######################################
        # get artifacts for the dependencies #
        ######################################
        subprocess.check_call(['mkdir', '-p', 'deps'])
        subprocess.check_call(['mkdir', '-p', 'deps/highs'])
        
        VERSION_BQPD = "1.0.0"
        VERSION_MUMPS = "5.8.0"
        VERSION_HSL = "4.0.4"
        VERSION_HIGHS = "1.11.0"
        
        artifact_version = get_machine() + '-' + get_system() + "-libgfortran5"
        
        # BQPD
        if not os.path.exists(current_directory + "/deps/lib/libbqpd.a"):
            bqpd_filename = "BQPD.v" + VERSION_BQPD + "." + artifact_version + ".tar.gz"
            subprocess.check_call(['wget', "-N", "https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-v" + VERSION_BQPD + "%2B0/" + bqpd_filename])
            subprocess.check_call(['tar', '-xzvf', bqpd_filename, '-C', 'deps'])
        
        # MUMPS
        if not os.path.exists(current_directory + "/deps/lib/libdmumps.a"):
            mumps_filename = "MUMPS_static.v" + VERSION_MUMPS + "." + artifact_version + ".tar.gz"
            subprocess.check_call(['wget', "-N", "https://github.com/amontoison/MUMPS_static_jll.jl/releases/download/MUMPS_static-v" + VERSION_MUMPS + "%2B0/" + mumps_filename])
            subprocess.check_call(['tar', '-xzvf', mumps_filename, '-C', 'deps'])
        
        # HSL
        if not os.path.exists(current_directory + "/deps/lib/libhsl.so"):
            hsl_filename = "HSL.v" + VERSION_HSL + "." + artifact_version + ".tar.gz"
            subprocess.check_call(['wget', "-N", "https://github.com/JuliaBinaryWrappers/HSL_jll.jl/releases/download/HSL-v" + VERSION_HSL + "%2B0/" + hsl_filename])
            subprocess.check_call(['tar', '-xzvf', hsl_filename, '-C', 'deps'])
        
        # HiGHS
        if not os.path.exists(current_directory + "/deps/highs/lib/cmake/highs"):
            highs_filename = "HiGHS_static.v" + VERSION_HIGHS + "." + artifact_version + ".tar.gz"
            subprocess.check_call(['wget', "-N", "https://github.com/amontoison/HiGHS_static_jll.jl/releases/download/HiGHS_static-v" + VERSION_HIGHS + "%2B0/" + highs_filename])
            subprocess.check_call(['tar', '-xzvf', highs_filename, '-C', 'deps/highs'])

        ##############################
        # compile the shared library #
        ##############################
        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            f"-DBQPD=" + current_directory + "/deps/lib/libbqpd.a",
            f"-DMUMPS_INCLUDE_DIR=" + current_directory + "/deps/include",
            f"-DMETIS_INCLUDE_DIR=" + current_directory + "/deps/include",
            f"-DMETIS_LIBRARY=" + current_directory + "/deps/lib/libmetis.a",
            f"-DMUMPS_LIBRARY=" + current_directory + "/deps/lib/libdmumps.a",
            f"-DMUMPS_COMMON_LIBRARY=" + current_directory + "/deps/lib/libmumps_common.a",
            f"-DMUMPS_PORD_LIBRARY=" + current_directory + "/deps/lib/libpord.a",
            f"-DMUMPS_MPISEQ_LIBRARY=" + current_directory + "/deps/lib/libmpiseq.a",
            f"-DBLAS_LIBRARIES=" + current_directory + "/deps/lib/libblas.a",
            f"-DLAPACK_LIBRARIES=" + current_directory + "/deps/lib/liblapack.a",
            #f"-DHSL=" + current_directory + "/deps/lib/libhsl.so",
            #f"-DHIGHS_DIR=" + current_directory + "/deps/highs/lib/cmake/highs",
        ]
        build_args = ["--target", "unopy"]
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # In this example, we pass in the version to C++. You might not need to.
        # cmake_args += [f"-DEXAMPLE_VERSION_INFO={self.distribution.get_version()}"]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator:
                try:
                    import ninja  # noqa: F401

                    cmake_args += ["-GNinja"]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )

setup(
    name="unopy",
    version="2.0.3",
    description="Uno: a next-gen Lagrange-Newton solver for nonlinearly constrained optimization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Charlie Vanaret",
    author_email="vanaret@zib.de",
    url="https://github.com/cvanaret/Uno",
    download_url="https://github.com/cvanaret/Uno/releases",
    project_urls={
        "Bug Tracker": "https://github.com/cvanaret/Uno/issues/",
        "Source Code": "https://github.com/cvanaret/Uno",
    },
    license="MIT",
    keywords="mathematics, optimization",
    ext_modules=[CMakeExtension("unopy")],
    cmdclass={
		"build_ext": CMakeBuild
	 },
    install_requires=["pybind11[global]"],
    python_requires=">=3.6",
    zip_safe=False,
)