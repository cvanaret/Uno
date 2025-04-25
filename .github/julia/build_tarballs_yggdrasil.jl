# Note that this script can accept some limited command-line arguments, run
# `julia build_tarballs.jl --help` to see a usage message.
using BinaryBuilder, Pkg

haskey(ENV, "UNO_RELEASE") || error("The environment variable UNO_RELEASE is not defined.")
haskey(ENV, "UNO_COMMIT") || error("The environment variable UNO_COMMIT is not defined.")
haskey(ENV, "UNO_URL") || error("The environment variable UNO_URL is not defined.")

name = "Uno"
version = VersionNumber(ENV["UNO_RELEASE"])

# Collection of sources required to complete build
sources = [
    GitSource(ENV["UNO_URL"], ENV["UNO_COMMIT"]),
]

# The remainder of the file is similar to the content found at the following link:
# https://raw.githubusercontent.com/JuliaPackaging/Yggdrasil/refs/heads/master/U/Uno/build_tarballs.jl
script = raw"""
cd $WORKSPACE/srcdir/Uno
mkdir -p build
cd build

if [[ "${target}" == *mingw* ]]; then
    LBT=blastrampoline-5
    HIGHS_DIR=${prefix}/lib
else
    LBT=blastrampoline
    HIGHS_DIR=${libdir}
fi

if [[ "${target}" == *apple* ]] || [[ "${target}" == *freebsd* ]]; then
    OMP=omp
else
    OMP=gomp
fi

cmake \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_PREFIX_PATH=${libdir} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Debug \
    -DAMPLSOLVER=${libdir}/libasl.${dlext} \
    -DHIGHS_DIR=${HIGHS_DIR} \
    -DBQPD=${prefix}/lib/libbqpd.a \
    -DHSL=${libdir}/libhsl.${dlext} \
    -DBLA_VENDOR="libblastrampoline" \
    -DMUMPS_INCLUDE_DIR=${includedir} \
    -DMETIS_INCLUDE_DIR=${includedir} \
    -DMUMPS_LIBRARY="${libdir}/libdmumps.${dlext}" \
    -DMUMPS_COMMON_LIBRARY="${libdir}/libmumps_common.${dlext}" \
    -DMUMPS_PORD_LIBRARY="${libdir}/libpord.${dlext}" \
    -DMUMPS_MPISEQ_LIBRARY="${libdir}/libmpiseq.${dlext}" \
    -DBLAS_LIBRARIES="${libdir}/lib${LBT}.${dlext}" \
    -DLAPACK_LIBRARIES="${libdir}/lib${LBT}.${dlext}" \
    -DBUILD_STATIC_LIBS=ON \
    -DBUILD_SHARED_LIBS=ON \
    ..

make uno_ampl -j${nproc}
make install

# Uno
install_license ${WORKSPACE}/srcdir/Uno/LICENSE

# BQPD
cp ${prefix}/share/licenses/BQPD/LICENSE ${WORKSPACE}/srcdir/Uno/LICENSE_BQPD
install_license ${WORKSPACE}/srcdir/Uno/LICENSE_BQPD
"""

platforms = supported_platforms()
filter!(p -> !(Sys.isfreebsd(p) && arch(p) == "aarch64"), platforms)
platforms = expand_cxxstring_abis(platforms)
platforms = expand_gfortran_versions(platforms)

products = [
    # This LibraryProduct may be useful once Uno provides a C API. We omit it for now.
    # LibraryProduct("libuno", :libuno),
    # We call this amplexe to match the convention of other JLL packages (like Ipopt_jll) that provide AMPL wrappers
    ExecutableProduct("uno_ampl", :amplexe),
    LibraryProduct("libuno", :libuno),
    FileProduct("lib/libuno.a", :libuno_a),
]

dependencies = [
    BuildDependency(PackageSpec(name="BQPD_jll", uuid="1325ac01-0a49-589f-8355-43321054aaab")),
    Dependency(PackageSpec(name="HiGHS_jll", uuid="8fd58aa0-07eb-5a78-9b36-339c94fd15ea"), compat="1.11.0"),
    Dependency(PackageSpec(name="HSL_jll", uuid="017b0a0e-03f4-516a-9b91-836bbd1904dd")),
    Dependency(PackageSpec(name="METIS_jll", uuid="d00139f3-1899-568f-a2f0-47f597d42d70")),
    Dependency(PackageSpec(name="ASL_jll", uuid="ae81ac8f-d209-56e5-92de-9978fef736f9"), compat="0.1.3"),
    Dependency(PackageSpec(name="MUMPS_seq_jll", uuid="d7ed1dd3-d0ae-5e8e-bfb4-87a502085b8d")),
    Dependency(PackageSpec(name="libblastrampoline_jll", uuid="8e850b90-86db-534c-a0d3-1478176c7d93"), compat="5.4.0"),
    # For OpenMP we use libomp from `LLVMOpenMP_jll` where we use LLVM as compiler (BSD systems),
    # and libgomp from `CompilerSupportLibraries_jll` everywhere else.
    Dependency(PackageSpec(name="CompilerSupportLibraries_jll", uuid="e66e0078-7015-5450-92f7-15fbd957f2ae"); platforms=filter(!Sys.isbsd, platforms)),
    Dependency(PackageSpec(name="LLVMOpenMP_jll", uuid="1d63c593-3942-5779-bab2-d838dc0a180e"); platforms=filter(Sys.isbsd, platforms)),

    # We need at least 3.29 (Ygg version), or 3.30 upstream version for LBT support,
    # so always pull the most recent CMake version.
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
    julia_compat = "1.9",
    preferred_gcc_version = v"10.2.0",
    clang_use_lld=false,
)