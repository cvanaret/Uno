## Installation instructions

### Packages and libraries

* download the AMPL solver library (ASL): http://www.netlib.org/ampl/solvers/

* download **optional** solvers:
    * BQPD (indefinite null-space QP solver): https://www.mcs.anl.gov/~leyffer/solvers.html
    * MA57 (sparse indefinite symmetric linear solver): http://www.hsl.rl.ac.uk/catalogue/ma57.html
    * LIBHSL (collection of libraries for sparse linear systems): https://licences.stfc.ac.uk/product/libhsl
    * MUMPS (sparse indefinite symmetric linear solver): https://mumps-solver.org/index.php?page=dwnld
    * HiGHS (linear programming solver): https://highs.dev

* to compile MUMPS in sequential mode, set the following variables at the end of your Makefile.inc:
```console
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
```

* install BLAS and LAPACK:
```console
sudo apt-get install libblas-dev liblapack-dev
```
* install cmake (and optionally ccmake, CMake curses interface):
```console
sudo apt-get install cmake cmake-curses-gui
```

### Compilation

1. Create a `build` directory in the main directory:
```console
mkdir build
```
2. Move to the build directory:
```console
cd build/
```
3. Execute cmake (you may provide the paths to the libraries ASL, BQPD and MA57):  
```console
cmake -DBQPD=path -DMA57=path -DAMPLSOLVER=path -DCMAKE_BUILD_TYPE=[Release|Debug] ..
```
4. **(or)** Use ccmake to provide the paths to the required and optional libraries:
```console
ccmake ..
```
5. Compile (in parallel: `n` being the number of threads, e.g. 6):
```console
make -jn
```

To compile the code with different configurations, simply create a `build` directory for each configuration and perform instructions 1 to 5.

### Unit tests

6. Install the GoogleTest suite:
```console
sudo apt-get install googletest
```
7. Perform steps 2 and 3 with the flag
```console
-DWITH_GTEST=ON
```
8. Run the test suite:
```console
./run_unotest
```

### Autocompletion

To benefit from autocompletion, install the file `uno_ampl-completion.bash`:
```console
sudo cp uno_ampl-completion.bash /etc/bash_completion.d/
```
and open a new terminal.

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
