## Installation instructions

### Packages and libraries

* download the AMPL solver library (ASL): http://www.netlib.org/ampl/solvers/

* **(optional)** download  solvers:
    * BQPD (null-space active set solver for nonconvex quadratic programming): get a precompiled binary for your architecture (https://github.com/leyffer/BQPD_jll.jl/releases) or get in touch with Sven Leyffer to apply for an academic license (https://www.mcs.anl.gov/~leyffer/solvers.html)
    * MA57 (sparse indefinite symmetric linear solver): http://www.hsl.rl.ac.uk/catalogue/ma57.html
    * LIBHSL (collection of libraries for sparse linear systems): https://licences.stfc.ac.uk/products/Software/HSL/LibHSL
    * MUMPS (sparse indefinite symmetric linear solver): https://mumps-solver.org/index.php?page=dwnld
    * HiGHS (linear programming solver): https://highs.dev

* to compile MUMPS in sequential mode, remove the flag `-fopenmp` at the end of your `Makefile.inc` and set the following variables:
```console
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
```

* **(optional)** install BLAS and LAPACK:
```console
sudo apt install libblas-dev liblapack-dev
```
* install cmake (and optionally ccmake, CMake curses interface):
```console
sudo apt install cmake cmake-curses-gui
```

### Compilation

1. Create a `build` directory in the main directory and move to it:
```console
mkdir build && cd build
```
2. Execute cmake:  
```console
cmake [options] ..
```
You can pass the following options:
- build type: `-DCMAKE_BUILD_TYPE=[Release|Debug]`
- build the Uno static library `uno_static`: `-DBUILD_STATIC_LIBS=[ON|OFF]`
- build the Uno shared library `uno_shared`: `-DBUILD_SHARED_LIBS=[ON|OFF]`
- enable LAPACK: `-DWITH_LAPACK=[ON|OFF]`
- path to the BQPD library: `-DBQPD=path_to_bqpd_lib`
- path to the MA27 library: `-DMA57=path_to_MA27_lib`
- path to the MA57 library: `-DMA57=path_to_MA57_lib`
- path to ASL library: `-DAMPLSOLVER=path_to_amplsolver_lib`
- path to HiGHS install directory: `-DHIGHS_DIR=path_to_highs_install_dir`
- path to HSL library: `-DHSL=path_to_hsl_lib`
- path to METIS library (`fakemetis` is built with MA57): `-DMETIS=path_to_metis_lib`
- path to MUMPS library: `-DMUMPS_LIBRARY=path_to_mumps_lib`
- path to MUMPS common library: `-DMUMPS_COMMON_LIBRARY=path_to_mumps_common_lib`
- path to MUMPS PORD library: `-DMUMPS_PORD_LIBRARY=path_to_mumps_pord_lib`
- path to MUMPS MPISEQ library: `-DMUMPS_MPISEQ_LIBRARY=path_to_mumps_mpiseq_lib`
- path to MUMPS include directory: `-DMUMPS_INCLUDE_DIR=path_to_mumps_include_dir`

3. **(or)** Use ccmake to provide the paths to the required and optional libraries:
```console
ccmake ..
```
4. Compile (in parallel: `n` being the number of threads, e.g. 6):
```console
make -jn
```
5. Compile the ASL executable:
```console
make uno_ampl -jn
```

To compile the code with different configurations, simply create a `build` directory for each configuration and perform instructions 1 to 5.

### Install

6. Install the built libraries and executable (`uno_static`, `uno_shared`, and `uno_ampl`):
```console
sudo make install
```

### Unit tests

7. Install the GoogleTest suite:
```console
sudo apt install googletest
```
8. Compile the test suite:
```console
make run_unotest -jn
```
9. Run the test suite:
```console
./run_unotest
```

### Precompiled libraries and executables

We provide precompiled Uno libraries and executables in the [releases tab](https://github.com/cvanaret/Uno/releases/latest/) for Linux, macOS (Intel & Silicon), and Windows.

On some platforms, the dynamic linker needs to know where to look for libraries at runtime.
You might need to set the following environment variables:

- `LD_LIBRARY_PATH` on Linux
- `DYLD_LIBRARY_PATH` or `DYLD_FALLBACK_LIBRARY_PATH` on macOS
- `PATH` on Windows

These variables should include the directory where you extracted the library files.
For all platforms, the environment variable `PATH` is needed to locate the binary `uno_ampl` / `uno_ampl.exe`.

**Example for Linux**:
```console
tar -xzf Uno.vX.Y.Z.linux.tar.gz
export LD_LIBRARY_PATH=/path/to/extracted/Uno/lib:$LD_LIBRARY_PATH
export PATH=/path/to/extracted/Uno/bin:$PATH
```

Note: The provided shared library `libhsl.so` / `libhsl.dylib` / `libhsl.dll` in the precompiled archive does not contain the HSL solvers like MA57 but can be replaced with the official version without the need to recompile anything.
