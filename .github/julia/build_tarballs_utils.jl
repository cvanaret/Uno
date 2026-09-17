# julia +1.7 --color=yes build_tarballs_utils.jl x86_64-linux-gnu-libgfortran5-cxx11,x86_64-apple-darwin-libgfortran5-cxx11,x86_64-w64-mingw32-libgfortran5-cxx11,aarch64-linux-gnu-libgfortran5-cxx11,aarch64-apple-darwin-libgfortran5-cxx11 --verbose --deploy="amontoison/UnoUtils_jll.jl"
#
# Note: the following variables is needed for cross-compilation on Mac platforms in the .bashrc
# export BINARYBUILDER_AUTOMATIC_APPLE="true"
using BinaryBuilder, Pkg

name = "UnoUtils"
version = v"2026.9.18"

# Collection of sources
sources = [
    # METIS v3.1.0
    GitSource("https://github.com/amontoison/METIS.git",
              "e827ffed17d56a4ac1add9cc33342c453a06c209"),
    FileSource("https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/refs/heads/master/M/METIS/METIS%405/bundled/patches/0001-mingw-w64-does-not-have-sys-resource-h.patch",
               "0ce0028dcb2205856aaf6f811a2fb80ec3d10290dd146cd3f060ba758de8d96f"),
    FileSource("https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/refs/heads/master/M/METIS/METIS%405/bundled/patches/0002-mingw-w64-do-not-use-reserved-double-underscored-names.patch",
               "5b49cc2f2d35f35c5946bc24d1a3ac408a09e6763b7a3250a6db5256dfea8293"),
    FileSource("https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/refs/heads/master/M/METIS/METIS%405/bundled/patches/0003-WIN32-Install-RUNTIME-to-bin.patch",
               "643ff86b8c587f718b0ba45d8a7f9d0f45484d32f2bbf7832f9350363d0998bd"),
    FileSource("https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/refs/heads/master/M/METIS/METIS%405/bundled/patches/0004-Fix-GKLIB_PATH-default-for-out-of-tree-builds.patch",
               "7d5977fc16d29bb0dc584b974e2270acdf34125e5f8a474588a4ec6ed57f1c9a"),
    FileSource("https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/refs/heads/master/M/METIS/METIS%405/bundled/patches/005-add-ifndefs.patch",
               "714dbe001a50882779e528f6170eca36f14a92df969dea8b2fbeef57026685b7"),    
    # BLAS / LAPACK v3.12.1
    GitSource("https://github.com/Reference-LAPACK/lapack.git",
              "6ec7f2bc4ecf4c4a93496aa2fa519575bc0e39ca"),
    # OpenBLAS v0.3.34 (release tarball + sha used by the known-good Yggdrasil OpenBLAS_jll build,
    # instead of an untested bare commit)
    ArchiveSource("https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.34/OpenBLAS-0.3.34.tar.gz",
                  "cd7e129868320cc2d033afa920e31202dfe0b8066a5b66661900ccc0f197dfed"),
    # MUMPS v5.9.1
    ArchiveSource("https://mumps-solver.org/MUMPS_5.9.1.tar.gz",
                  "659c9b57646b5a003ac618baa1faf9dd2044e46c732b3daaccbc7158003e1b46"),
    # HiGHS v1.15.1
    GitSource("https://github.com/ERGO-Code/HiGHS.git",
              "04024d701f79feb8e2f18bc3df0dffc04ef05088"),
    # SPRAL v2025.9.18
    GitSource("https://github.com/ralna/spral.git",
              "80bc843ac3847d4a783a0e11213715a70175aee6"),
    # Hwloc v2.13.0
    # ArchiveSource("https://download.open-mpi.org/release/hwloc/v2.13/hwloc-2.13.0.tar.bz2",
    #               "52e936afb6ebd80f171f763fcf14f7b1f5ce98b125af5dd2f328b873b1fd0dab"),
    # Package compiler for Windows
    ArchiveSource("https://github.com/JuliaLang/PackageCompiler.jl/releases/download/v1.0.0/x86_64-8.1.0-release-posix-seh-rt_v6-rev0.tar.gz",
                  "fe3f401bc936fbe6af940b26c5e0f266f762a3416f979c706e599b24082dc5c7"),
]

# Bash recipe for building across all platforms
script = raw"""
# Remove system CMake to use the jll version
apk del cmake

# Update Ninja
cp ${host_prefix}/bin/ninja /usr/bin/ninja

## ----- Compile METIS -----
cd $WORKSPACE/srcdir/METIS
if [ $target = "x86_64-w64-mingw32" ] || [ $target = "i686-w64-mingw32" ]; then
    atomic_patch -p1 $WORKSPACE/srcdir/0001-mingw-w64-does-not-have-sys-resource-h.patch
    atomic_patch -p1 $WORKSPACE/srcdir/0002-mingw-w64-do-not-use-reserved-double-underscored-names.patch
    atomic_patch -p1 $WORKSPACE/srcdir/0003-WIN32-Install-RUNTIME-to-bin.patch
    atomic_patch -p1 $WORKSPACE/srcdir/0004-Fix-GKLIB_PATH-default-for-out-of-tree-builds.patch
fi
atomic_patch -p1 $WORKSPACE/srcdir/005-add-ifndefs.patch

# fix CMake version
sed -i 's/VERSION 2.8/VERSION 3.5/' CMakeLists.txt

mkdir -p build
cd build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DCMAKE_TOOLCHAIN_FILE="${CMAKE_TARGET_TOOLCHAIN}" \
    -DCMAKE_VERBOSE_MAKEFILE=1 \
    -DGKLIB_PATH=$WORKSPACE/srcdir/METIS/GKlib \
    -DSHARED=0
make -j${nproc}
make install

## ----- Compile BLAS / LAPACK -----
cd $WORKSPACE/srcdir/lapack
mkdir build
cd build
cmake .. \
 -DCBLAS=ON \
 -DLAPACKE=OFF \
 -DCMAKE_INSTALL_PREFIX=${prefix} \
 -DCMAKE_TOOLCHAIN_FILE="${CMAKE_TARGET_TOOLCHAIN}" \
 -DCMAKE_BUILD_TYPE=Release \
 -DBUILD_SHARED_LIBS=OFF \
 -DBUILD_INDEX64_EXT_API=OFF \
 -DTEST_FORTRAN_COMPILER=OFF
make -j$(nproc)
make install

## ----- Compile OpenBLAS -----
cd $WORKSPACE/srcdir/OpenBLAS*

# Fix the issue with the symbol __imp__cprintf on Windows.
if [[ ${target} == *mingw* ]]; then
    for f in driver/others/memory.c driver/others/xerbla.c; do
        sed -i '/^#define[[:space:]]\+printf[[:space:]]\+_cprintf[[:space:]]*$/d' "${f}"
        if grep -q '_cprintf' "${f}"; then
            echo "ERROR: _cprintf still referenced in ${f}"
            exit 1
        fi
    done
fi

# We always want threading
flags=(USE_THREAD=1 GEMM_MULTITHREADING_THRESHOLD=400 NO_AFFINITY=1)

# We are cross-compiling
flags+=(CROSS=1 "CROSS_SUFFIX=${target}-")

# We need to use our basic objconv, not a prefixed one
flags+=(OBJCONV=objconv)

# Static library only, with the 32-bit integer (LP64) interface
flags+=(NO_SHARED=1 INTERFACE64=0 LIBPREFIX=libopenblas)

# Word size and maximum thread count
if [[ ${nbits} == 32 ]]; then
    flags+=(BINARY=32 NUM_THREADS=8)
else
    flags+=(NUM_THREADS=32)
fi
if [[ ${target} == x86_64-* ]]; then
    flags+=(BINARY=64)
fi

# Runtime kernel dispatch. We ship these binaries to unknown hardware,
# so we embed every kernel set and let OpenBLAS pick at runtime.
if [[ ${proc_family} == intel ]]; then
    flags+=(DYNAMIC_ARCH=1 TARGET=GENERIC)
elif [[ ${target} == aarch64-* ]]; then
    flags+=(TARGET=ARMV8 DYNAMIC_ARCH=1)
elif [[ ${target} == arm-* ]]; then
    flags+=(TARGET=ARMV7)
elif [[ ${target} == powerpc64le-* ]]; then
    flags+=(TARGET=POWER8 DYNAMIC_ARCH=1)
elif [[ ${target} == riscv64-* ]]; then
    flags+=(TARGET=RISCV64_GENERIC DYNAMIC_ARCH=1)
fi

# SME is supported neither by pre-M4 hardware nor by our Darwin toolchains
if [[ ${target} == aarch64-*-darwin* ]]; then
    export NO_SME=1
fi

if [[ ${target} == x86_64-w64-mingw32 ]]; then
    flags+=("CFLAGS=${CFLAGS} -fno-asynchronous-unwind-tables")
fi

# Choose our make parallelism. The Makefile would otherwise override our choice.
flags+=(-j${nproc})
export MAKE_NB_JOBS=0

echo "OpenBLAS build flags: ${flags[@]}"
make "${flags[@]}"
make "${flags[@]}" "PREFIX=$prefix" install

cd ${prefix}/lib
if [[ ! -f libopenblas.a || -L libopenblas.a ]]; then
    versioned_a=$(ls libopenblas*-r*.a 2>/dev/null | head -1)
    if [[ -n "${versioned_a}" ]]; then
        rm -f libopenblas.a
        mv "${versioned_a}" libopenblas.a
    fi
fi
rm -f libopenblas*-r*.a
ls -la ${prefix}/lib

## ----- Compile MUMPS -----
cd $WORKSPACE/srcdir/MUMPS*

makefile="Makefile.G95.SEQ"
cp Make.inc/${makefile} Makefile.inc

# Add `-fallow-argument-mismatch` if supported
: >empty.f
FFLAGS=()
if gfortran -c -fallow-argument-mismatch empty.f >/dev/null 2>&1; then
    FFLAGS+=("-fallow-argument-mismatch")
fi
rm -f empty.*

make_args+=(OPTF="-fPIC -O3"
            OPTL="-fPIC -O3"
            OPTC="-fPIC -O3"
            CDEFS=-DAdd_
            LMETISDIR=${libdir}
            IMETIS=-I${includedir}
            LMETIS="-L${libdir} -lmetis"
            ORDERINGSF="-Dpord -Dmetis"
            CC="$CC ${CFLAGS[@]}"
            FC="gfortran ${FFLAGS[@]}"
            FL="gfortran"
            RANLIB="echo"
            LIBBLAS="-L${libdir} -lopenblas"
            LAPACK="-L${libdir} -lopenblas")

make -j${nproc} d "${make_args[@]}"
cp include/*.h ${includedir}
cp lib/*.a ${prefix}/lib

## ----- Compile Hwloc -----
# cd $WORKSPACE/srcdir/hwloc-*
# if [[ "${target}" == *-apple-darwin* ]]; then
#     ./configure --prefix=${prefix} --build=${MACHTYPE} --host=${target} --disable-static --enable-shared
# else
#     CFLAGS="${CFLAGS} -fPIC" ./configure --prefix=${prefix} --build=${MACHTYPE} --host=${target} --enable-static --disable-shared
# fi
# make -j${nproc}
# make install

## ----- Compile SPRAL -----
cd $WORKSPACE/srcdir/spral

meson setup builddir --cross-file="${MESON_TARGET_TOOLCHAIN}" \
                     --prefix=$prefix \
                     -Ddefault_library=static \
                     -Dlibhwloc= \
                     -Dmodules=false \
                     -Dopenmp=false \
                     -Dlibblas=openblas \
                     -Dliblapack=openblas \
                     -Dbinaries=false \
                     -Dtests=false \
                     -Dexamples=false

meson compile -C builddir
meson install -C builddir

## ----- Compile HiGHS -----
cd $WORKSPACE/srcdir/HiGHS

# Patch v1.15.1 (see https://github.com/JuliaPackaging/Yggdrasil/tree/master/H/HiGHS/bundled/patches)
# fix-cli11.patch
sed -i 's/(*opt)/opt->count() > 0/' extern/cli11/CLI11.hpp
# fix-destroy.patch
sed -i 's/Highs::resetGlobalScheduler(true);//' highs/interfaces/highs_c_api.cpp

# On macOS, HiGHS hard-codes `-framework Accelerate` for HiPO in
# cmake/FindHipoDeps.cmake and consults neither BLA_VENDOR nor BLAS_LIBRARIES.
# Neutralise the three APPLE branches (highs_configure_blas_target,
# highs_configure_blas_metadata, highs_link_blas) so HiPO links the OpenBLAS we
# ship: a single BLAS in the binary rather than Accelerate alongside OpenBLAS.
# HIPO_USES_APPLE_BLAS is licensing metadata only, no C++ source depends on it.
if [[ "${target}" == *apple* ]]; then
    sed -i 's/^    if(APPLE)$/    if(FALSE) # Uno: link our OpenBLAS, not Accelerate/' cmake/FindHipoDeps.cmake
    # Fail loudly instead of silently falling back to Accelerate
    grep -q 'if(FALSE) # Uno' cmake/FindHipoDeps.cmake || \
        { echo "ERROR: the HiGHS Accelerate patch did not apply"; exit 1; }
fi

mkdir build
cd build
cmake .. \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF \
    -DZLIB=OFF \
    -DHIPO=ON \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_TESTING=OFF \
    -DBUILD_CXX_EXE=OFF \
    -DBLA_VENDOR=Generic \
    -DBLAS_LIBRARIES=${prefix}/lib/libblas.a \
    -DBUILD_SHARED_EXTRAS_LIB=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON

if [[ "${target}" == *-linux* ]]; then
        make -j ${nproc}
else
    if [[ "${target}" == *-mingw* ]]; then
        cmake --build . --config Release
    else
        cmake --build . --config Release --parallel
    fi
fi
make install

## ----- Second flavour of highs_extras, this one tied to OpenBLAS -----
cd $WORKSPACE/srcdir/HiGHS
mkdir -p build_openblas
cd build_openblas
cmake .. \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF \
    -DZLIB=OFF \
    -DHIPO=ON \
    -DBUILD_EXAMPLES=OFF \
    -DBUILD_TESTING=OFF \
    -DBUILD_CXX_EXE=OFF \
    -DBLA_VENDOR=OpenBLAS \
    -DBLAS_LIBRARIES=${prefix}/lib/libopenblas.a \
    -DOPENBLAS_LIB=${prefix}/lib/libopenblas.a \
    -DOPENBLAS_INCLUDE_DIR=${includedir} \
    -DBUILD_SHARED_EXTRAS_LIB=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON

cmake --build . --target highs_extras --parallel ${nproc}

openblas_extras=$(find . -name 'libhighs_extras.a' | head -1)
if [[ -z "${openblas_extras}" ]]; then
    echo "ERROR: libhighs_extras.a not found in the OpenBLAS HiGHS build"
    exit 1
fi
cp "${openblas_extras}" ${prefix}/lib/libhighs_extras_openblas.a

if [[ "${target}" == *-mingw* ]]; then
    cp $WORKSPACE/srcdir/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/libstdc++.a ${prefix}/lib/libstdc++.a
    cp $WORKSPACE/srcdir/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/libgfortran.a ${prefix}/lib/libgfortran.a
    cp $WORKSPACE/srcdir/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/libquadmath.a ${prefix}/lib/libquadmath.a
    cp $WORKSPACE/srcdir/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/libgomp.a ${prefix}/lib/libgomp.a
    cp $WORKSPACE/srcdir/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/libgcc.a ${prefix}/lib/libgcc.a
    cp $WORKSPACE/srcdir/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/libgcc_eh.a ${prefix}/lib/libgcc_eh.a
fi

# Clean
rm -r ${prefix}/bin
rm ${prefix}/lib/libsmumps.a

# Compile some shared libraries for Windows
if [ $target = "x86_64-w64-mingw32" ] || [ $target = "i686-w64-mingw32" ]; then
    ## METIS
    cd $WORKSPACE/srcdir/METIS
    mkdir -p build_shared
    cd build_shared
    cmake .. \
        -DCMAKE_INSTALL_PREFIX=$prefix \
        -DCMAKE_TOOLCHAIN_FILE="${CMAKE_TARGET_TOOLCHAIN}" \
        -DCMAKE_VERBOSE_MAKEFILE=1 \
        -DGKLIB_PATH=$WORKSPACE/srcdir/METIS/GKlib \
        -DSHARED=1
    make -j${nproc}
    make install
fi
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = supported_platforms()
platforms = expand_gfortran_versions(platforms)

# The products that we will ensure are always built
products = [
    FileProduct("lib/libmetis.a", :libmetis_a),
    FileProduct("lib/libblas.a", :libblas_a),
    FileProduct("lib/libcblas.a", :libcblas_a),
    FileProduct("lib/liblapack.a", :liblapack_a),
    FileProduct("lib/libopenblas.a", :libopenblas_a),
    FileProduct("lib/libpord.a", :libpord_a),
    FileProduct("lib/libmpiseq.a", :libmpiseq_a),
    FileProduct("lib/libmumps_common.a", :libmumps_common_a),
    FileProduct("lib/libdmumps.a", :libdmumps_a),
    FileProduct("lib/libhighs.a", :libhighs_a),
    FileProduct("lib/libhighs_extras.a", :libhighs_extras_a),
    FileProduct("lib/libhighs_extras_openblas.a", :libhighs_extras_openblas_a),
    FileProduct("lib/libspral.a", :libspral_a),
    # FileProduct("lib/libhwloc.a", :libhwloc_a),
]

# Dependencies that must be installed before this package can be built
dependencies = [
    HostBuildDependency(PackageSpec(name="Ninja_jll", uuid="76642167-d241-5cee-8c94-7a494e8cb7b7")),
    HostBuildDependency(PackageSpec(name="CMake_jll", uuid="3f4e10e2-61f2-5801-8945-23b9d642d0e6")),
    Dependency(PackageSpec(name="LLVMOpenMP_jll", uuid="1d63c593-3942-5779-bab2-d838dc0a180e"); platforms=filter(Sys.isbsd, platforms)),
]

build_tarballs(
    ARGS,
    name,
    version,
    sources,
    script,
    platforms,
    products,
    dependencies;
    julia_compat = "1.6",
    preferred_gcc_version = v"12",
    preferred_llvm_version = v"18.1.7",
    clang_use_lld=false,
    lock_microarchitecture=false,
)
