# julia +1.7 --color=yes build_tarballs_utils.jl x86_64-linux-gnu-libgfortran5,x86_64-apple-darwin-libgfortran5,x86_64-w64-mingw32-libgfortran5,aarch64-linux-gnu-libgfortran5,aarch64-apple-darwin-libgfortran5 --verbose --deploy="amontoison/UnoUtils_jll.jl"
using BinaryBuilder, Pkg

name = "UnoUtils"
version = v"2026.2.7"

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
    # OpenBLAS v0.3.31
    # ArchiveSource("https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.31/OpenBLAS-0.3.31.tar.gz",
    #               "6dd2a63ac9d32643b7cc636eab57bf4e57d0ed1fff926dfbc5d3d97f2d2be3a6"),
    # MUMPS v5.8.2
    ArchiveSource("https://mumps-solver.org/MUMPS_5.8.2.tar.gz",
                  "eb515aa688e6dbab414bb6e889ff4c8b23f1691a843c68da5230a33ac4db7039"),
    # HiGHS v1.13.0
    GitSource("https://github.com/ERGO-Code/HiGHS.git",
              "1bce6d5c801398dab6d2e6f98ac8935f3d4eec9c"),
]

# Bash recipe for building across all platforms
script = raw"""
# Remove system CMake to use the jll version
apk del cmake

## ----- Compile METIS -----
cd $WORKSPACE/srcdir/METIS
if [ $target = "x86_64-w64-mingw32" ] || [ $target = "i686-w64-mingw32" ]; then
    atomic_patch -p1 $WORKSPACE/srcdir/0001-mingw-w64-does-not-have-sys-resource-h.patch
    atomic_patch -p1 $WORKSPACE/srcdir/0002-mingw-w64-do-not-use-reserved-double-underscored-names.patch
    atomic_patch -p1 $WORKSPACE/srcdir/0003-WIN32-Install-RUNTIME-to-bin.patch
    atomic_patch -p1 $WORKSPACE/srcdir/0004-Fix-GKLIB_PATH-default-for-out-of-tree-builds.patch
fi
atomic_patch -p1 $WORKSPACE/srcdir/005-add-ifndefs.patch
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
# cd ${WORKSPACE}/srcdir/OpenBLAS*/

# # We always want threading
# flags=(USE_THREAD=1 GEMM_MULTITHREADING_THRESHOLD=400 NO_AFFINITY=1)
# if [[ "${CONSISTENT_FPCSR}" == "true" ]]; then
#     flags+=(CONSISTENT_FPCSR=1)
# fi

# # We are cross-compiling
# flags+=(CROSS=1 PREFIX=/ "CROSS_SUFFIX=${target}-")

# # We need to use our basic objconv, not a prefixed one:
# flags+=(OBJCONV=objconv)

# LIBPREFIX=libopenblas
# flags+=("LIBPREFIX=${LIBPREFIX}")

# # If we're building for x86_64 Windows gcc7+, we need to disable usage of
# # certain AVX-512 registers (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=65782)
# if [[ ${target} == x86_64-w64-mingw32 ]] && [[ $(gcc --version | head -1 | awk '{ print $3 }') =~ (7|8).* ]]; then
#     CFLAGS="${CFLAGS} -fno-asynchronous-unwind-tables"
# fi

# # Because we use this OpenBLAS within Julia, and often want to bundle our
# # libgfortran and other friends alongside, we need an RPATH of '$ORIGIN',
# # so set it here.
# if [[ ${target} == *linux* ]] || [[ ${target} == *freebsd* ]]; then
#     export LDFLAGS="${LDFLAGS} '-Wl,-rpath,\$\$ORIGIN' -Wl,-z,origin"
# elif [[ ${target} == *apple* ]]; then
#     export LDFLAGS="${LDFLAGS} -Wl,-rpath,@loader_path/"
# fi

# # Choose our make parallelism.
# flags+=(-j${nproc})

# # The Makefile will otherwise override our choice
# export MAKE_NB_JOBS=0

# # Build the actual library
# make "${flags[@]}"

# # Install the library
# make "${flags[@]}" "PREFIX=$prefix" install

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

make_args+=(OPTF="-O3"
            OPTL="-O3"
            OPTC="-O3"
            CDEFS=-DAdd_
            LMETISDIR=${libdir}
            IMETIS=-I${includedir}
            LMETIS="-L${libdir} -lmetis"
            ORDERINGSF="-Dpord -Dmetis"
            CC="$CC ${CFLAGS[@]}"
            FC="gfortran ${FFLAGS[@]}"
            FL="gfortran"
            RANLIB="echo"
            LIBBLAS="-L${libdir} -lblas"
            LAPACK="-L${libdir} -llapack")

make -j${nproc} d "${make_args[@]}"
cp include/*.h ${includedir}
cp lib/*.a ${prefix}/lib

## ----- Compile HiGHS -----
cd $WORKSPACE/srcdir/HiGHS
mkdir build
cd build
if [[ "${target}" == *-mingw* ]]; then
    LIBGFORTRAN=libgfortran-5
else
    LIBGFORTRAN=libgfortran
fi
cmake .. \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=OFF \
    -DHIPO=ON \
    -DBUILD_TESTING=OFF \
    -DUILD_CXX_EXE=OFF \
    -DBLAS_LIBRARIES="${prefix}/lib/libcblas.a;${prefix}/lib/libblas.a;${libdir}/${LIBGFORTRAN}.${dlext}"

if [[ "${target}" == *-linux-* ]]; then
        make -j ${nproc}
else
    if [[ "${target}" == *-mingw* ]]; then
        cmake --build . --config Release
    else
        cmake --build . --config Release --parallel
    fi
fi
make install
"""

# These are the platforms we will build for by default, unless further
# platforms are passed in on the command line
platforms = supported_platforms()
platforms = expand_gfortran_versions(platforms)

# The products that we will ensure are always built
products = [
    FileProduct("lib/libbqpd.a", :libbqpd_a),
    FileProduct("lib/libmetis.a", :libmetis_a),
    FileProduct("lib/libblas.a", :libblas_a),
    FileProduct("lib/libcblas.a", :libcblas_a),
    FileProduct("lib/liblapack.a", :liblapack_a),
    # FileProduct("lib/libopenblas.a", :libopenblas_a),
    FileProduct("lib/libpord.a", :libpord_a),
    FileProduct("lib/libmpiseq.a", :libmpiseq_a),
    FileProduct("lib/libmumps_common.a", :libmumps_common_a),
    FileProduct("lib/libdmumps.a", :libdmumps_a),
    FileProduct("lib/libhighs.a", :libhighs_a),
]

# Dependencies that must be installed before this package can be built
dependencies = [
    BuildDependency(PackageSpec(name="BQPD_jll", uuid="1325ac01-0a49-589f-8355-43321054aaab")),
    HostBuildDependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae")),
    HostBuildDependency(PackageSpec(name="CMake_jll", uuid="3f4e10e2-61f2-5801-8945-23b9d642d0e6")),
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
    preferred_gcc_version = v"10.2.0",
    clang_use_lld=false,
)
