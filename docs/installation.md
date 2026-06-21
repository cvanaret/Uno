# Installation guide

**Contents**

- [Precompiled libraries and executables](#precompiled-libraries-and-executables)
- [Dependencies](#dependencies)
- [Compiling and installing Uno from source](#compiling-and-installing-uno-from-source)
- [CMake options](#cmake-options)

## Precompiled libraries and executables

We provide precompiled Uno libraries and executables in the [releases tab](https://github.com/cvanaret/Uno/releases/latest/) for Linux (x64 and aarch64), macOS (x64 and aarch64), and Windows (x64).

> [!NOTE]
> The shared library `libhsl.so` / `libhsl.dylib` / `libhsl.dll` provided in the precompiled archive is a dummy version that does not include the official HSL linear solvers such as `MA27` or `MA57`.
> However, it can be safely replaced with the official precompiled [libHSL](https://licences.stfc.ac.uk/products/Software/HSL/LibHSL) library without the need to recompile anything.
> The routine symbols are identical, allowing seamless hot-swapping of the library.

On some platforms, the dynamic linker needs to know where to look for libraries at runtime.
You might need to set the following environment variables:

- `LD_LIBRARY_PATH` on Linux
- `DYLD_LIBRARY_PATH` or `DYLD_FALLBACK_LIBRARY_PATH` on macOS
- `PATH` on Windows

These variables should include the directories where you extracted the library files.
For all platforms, the environment variable `PATH` is needed to locate the binary `uno_ampl` / `uno_ampl.exe`:

=== "Linux"

    ```console
    tar -xzf Uno.vX.Y.Z.linux.tar.gz
    export LD_LIBRARY_PATH=/path/to/extracted/Uno/lib:/path/to/extracted/Uno/deps:$LD_LIBRARY_PATH
    export PATH=/path/to/extracted/Uno/bin:$PATH
    ```

=== "macOS"

    ```console
    tar -xzf Uno.vX.Y.Z.macos.tar.gz
    export DYLD_LIBRARY_PATH=/path/to/extracted/Uno/lib:/path/to/extracted/Uno/deps:$DYLD_LIBRARY_PATH
    export PATH=/path/to/extracted/Uno/bin:$PATH
    ```

=== "Windows (PowerShell)"

    ```console
    tar -xzf Uno.vX.Y.Z.windows.zip
    $env:PATH="C:\path\to\extracted\Uno\bin;C:\path\to\extracted\Uno\lib;C:\path\to\extracted\Uno\deps;$env:PATH"
    ```

=== "Windows (Command Prompt)"

    ```console
    tar -xzf Uno.vX.Y.Z.windows.zip
    set PATH=C:\path\to\extracted\Uno\bin;C:\path\to\extracted\Uno\lib;C:\path\to\extracted\Uno\deps;%PATH%
    ```

## Dependencies

* install cmake, BLAS and LAPACK:
```console
sudo apt update
sudo apt install cmake
sudo apt install libblas-dev liblapack-dev
```

* **(optional)** download the [AMPL solver library](http://www.netlib.org/ampl/solvers/) (ASL)

* **(optional)** download  solvers:
    * BQPD (null-space active set solver for nonconvex quadratic programming): get a [precompiled binary](https://github.com/leyffer/BQPD_jll.jl/releases) for your architecture or [get in touch with Sven Leyffer](https://www.mcs.anl.gov/~leyffer/solvers.html) to apply for an academic license
    * [MA27](https://www.hsl.rl.ac.uk/download/MA27/1.0.0) (sparse indefinite symmetric linear solver)
    * [MA57](http://www.hsl.rl.ac.uk/catalogue/ma57.html) (sparse indefinite symmetric linear solver)
    * [LIBHSL](https://licences.stfc.ac.uk/products/Software/HSL/LibHSL) (collection of solvers for sparse linear systems)
    * [MUMPS](https://mumps-solver.org/index.php?page=dwnld) (sparse indefinite symmetric linear solver)
    * [SSIDS](https://github.com/ralna/spral) (sparse indefinite symmetric linear solver) from the SPRAL library
    * [HiGHS](https://highs.dev) (linear programming and convex quadratic programming solver; it can also ship the HiPO sparse symmetric indefinite linear solver, exposed in Uno as the `HIPO` linear solver)

* the `HIPO` linear solver requires a HiGHS built with HiPO support (HiGHS option `-DHIPO=ON`, off by default). Since HiGHS does not install the FactorHighs C API header (and there is no reliable way to detect HiPO from the installed HiGHS), enable it explicitly when configuring Uno with `-DHIGHS_WITH_HIPO=ON` (HiGHS must be found). Uno vendors the FactorHighs C API header and resolves the symbols from the linked `libhighs`

* to compile MUMPS in sequential mode, remove the flag `-fopenmp` at the end of your `Makefile.inc` and set the following variables:
```console
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
```

* compile SSIDS without OpenMP with the Meson flag `-Dopenmp=false`

* you may experience a short lag at startup (about 1/4s) when running Uno with SSIDS. This is due to `hwloc` (hardware locality), a tool that aims at discovering hardware resources in parallel architectures. To precompute the required topology, run the following commands before running Uno:
```
lstopo --of xml ~/.config/hwloc-topology.xml
echo 'export HWLOC_XMLFILE=$HOME/.config/hwloc-topology.xml' >> ~/.bashrc
```

## Compiling and installing Uno from source

The sequence of commands to configure and build is as follows (assuming the build directory is `build`):
```console
cmake -S . -B build [options]
cmake --build build --parallel
```
See the list of CMake options [here](#cmake-options).

To install the built libraries and headers:
```console
cmake --install build
```

To configure, compile and run the test suite:
```console
cmake -S . -B build -DENABLE_TESTS=ON
cmake --build build --target run_unotest --parallel
ctest --test-dir build
```

## CMake options

You can pass the following options as `-DOPTION=value`:

| Option                 | Description                                                                                                | Possible values             |
|:-----------------------|:-----------------------------------------------------------------------------------------------------------|:----------------------------|
| `CMAKE_BUILD_TYPE`     | build type                                                                                                 | `Release`, `Debug`          |
| `ENABLE_TESTS`         | enable the unit tests                                                                                      | `ON`, `OFF`                 |
| `BUILD_STATIC_LIBS`    | build the Uno static library                                                                               | `ON`, `OFF`                 |
| `BUILD_SHARED_LIBS`    | build the Uno shared library                                                                               | `ON`, `OFF`                 |
| `LAPACK_LIBRARIES`     | path(s) to the LAPACK library, separated by `;`                                                            | `paths`                     |
| `BLAS_LIBRARIES`       | path(s) to the BLAS library, separated by `;`                                                              | `paths`                     |
| `AMPLSOLVER`           | path(s) to the ASL library                                                                                 | `path_to_libamplsolver`     |
| `BQPD`                 | path to the BQPD library                                                                                   | `path_to_libbqpd`           |
| `MA27`                 | path to the MA27 library                                                                                   | `path_to_libma27`           |
| `MA57`                 | path to the MA57 library                                                                                   | `path_to_libma57`           |
| `HSL`                  | path to the HSL library                                                                                    | `path_to_libhsl`            |
| `HIGHS`                | path to the HiGHS libraries (typically `libhighs` and `libhighs_extras`)                                   | `path_to_libhighs`          |
| `HIGHS_INCLUDE_DIR`    | path to HiGHS include directory                                                                            | `path_to_highs_include_dir` |
| `METIS`                | path to the METIS library                                                                                  | `path_to_libmetis`          |
| `MUMPS_LIBRARY`        | path to the MUMPS library                                                                                  | `path_to_libdmumps`         |
| `MUMPS_COMMON_LIBRARY` | path to the MUMPS common library                                                                           | `path_to_libmumps_common`   |
| `MUMPS_PORD_LIBRARY`   | path to the MUMPS PORD library                                                                             | `path_to_libpord`           |
| `MUMPS_MPISEQ_LIBRARY` | path to the MUMPS MPISEQ library                                                                           | `path_to_libmpiseq`         |
| `MUMPS_INCLUDE_DIR`    | path to MUMPS include directory                                                                            | `path_to_mumps_include_dir` |
| `SPRAL`                | path to SPRAL library </br> (requires `libhwloc` passed to `AUXILIARY_LIBRARIES`)                          | `path_to_libspral`          |
| `SPRAL_INCLUDE_DIR`    | path to SPRAL include directory                                                                            | `path_to_spral_include_dir` |
| `AUXILIARY_LIBRARIES`  | path(s) to additional libraries to link against, separated by `;` </br> (e.g., `libhwloc` and `libstdc++`) | `paths`                     |
