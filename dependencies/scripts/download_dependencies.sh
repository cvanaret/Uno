#!/usr/bin/env bash
set -e

# detect OS
OS_NAME="$(uname -s)"
case "$OS_NAME" in
    Linux*)
      if [[ "$CIBW_BUILD" == *musllinux* ]] || ldd --version 2>&1 | grep -qi musl; then
        OS="linux-musl"
      else
        OS="linux-gnu"
      fi
      ;;
    Darwin*) OS="apple-darwin";;
    Windows*) OS="w64-mingw32";;
    MINGW64_NT*) OS="w64-mingw32";;
    MSYS_NT*) OS="w64-mingw32";;
    *) echo "Unsupported OS: $OS_NAME"; exit 1;;
esac
# detect architecture
ARCH_NAME="$(uname -m)"
case "$ARCH_NAME" in
    x86_64|amd64) ARCH="x86_64";;
    arm64|aarch64) ARCH="aarch64";;
    *) echo "Unknown architecture '$ARCH_NAME'."; exit 1;;
esac

# change directory
cd dependencies

# download BQPD
VERSION="v1.0.0"
REPO="https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-${VERSION}%2B0"
ASSET_NAME="BQPD.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
curl -L -o BQPD.tar.gz "$ASSET_URL"
tar -xzf BQPD.tar.gz
pwd

# download UnoUtils: MUMPS (+ METIS, BLAS and LAPACK)
VERSION="2026.7.2"
REPO="https://github.com/amontoison/UnoUtils_jll.jl/releases/download/UnoUtils-v${VERSION}%2B0"
ASSET_NAME="UnoUtils.v${VERSION}.${ARCH}-${OS}-libgfortran5-cxx11.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
curl -L -o UnoUtils.tar.gz "$ASSET_URL"
tar -xzf UnoUtils.tar.gz
pwd

# on Windows, compile HiGHS on the fly and replace that of UnoUtils
if [[ "$OS" == "w64-mingw32" ]]; then
	# delete UnoUtils' HiGHS
	rm -rf lib/libhighs* include/highs include/Highs*.h bin/highs*

	# fetch HiGHS source
	BUILD_ROOT="$(mktemp -d)"
	mkdir "${BUILD_ROOT}/HiGHS"
	VERSION="v1.15.0"
	ASSET_URL="https://github.com/ERGO-Code/HiGHS/archive/refs/tags/${VERSION}.tar.gz"
	echo "Downloading: ${ASSET_URL}"
	curl -fL -o "${BUILD_ROOT}/HiGHS-src.tar.gz" "$ASSET_URL"
	tar -xzf "${BUILD_ROOT}/HiGHS-src.tar.gz" -C "${BUILD_ROOT}/HiGHS" --strip-components=1

	# build
	cmake -S "${BUILD_ROOT}/HiGHS" -B "${BUILD_ROOT}/build" \
		-G "MinGW Makefiles" \
		-DCMAKE_INSTALL_PREFIX="${BUILD_ROOT}/install" \
		-DCMAKE_BUILD_TYPE=Release \
		-DBUILD_SHARED_LIBS=OFF \
		-DZLIB=OFF \
		-DHIPO=ON \
		-DBUILD_EXAMPLES=OFF \
		-DBUILD_TESTING=OFF \
		-DBUILD_CXX_EXE=OFF \
		-DBLAS_LIBRARIES="${PWD}/lib/libblas.a" \
		-DBUILD_SHARED_EXTRAS_LIB=OFF \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON

	cmake --build "${BUILD_ROOT}/build" --config Release --parallel
	cmake --install "${BUILD_ROOT}/build" --config Release

	cp -a "${BUILD_ROOT}/install/lib/."     lib
	cp -a "${BUILD_ROOT}/install/include/." include
	rm -rf "${BUILD_ROOT}"

	CXX="${CXX:-g++}"
	LIBSTDCXX="$("$CXX" -print-file-name=libstdc++.a)"
	if [[ "$LIBSTDCXX" == "libstdc++.a" || ! -f "$LIBSTDCXX" ]]; then
	  echo "Could not locate libstdc++.a via $CXX"; exit 1
	fi
	echo "Using $LIBSTDCXX ($($CXX -dumpversion))"
	cp "$LIBSTDCXX" "lib/libstdc++.a"
fi

# delete unwanted directories
# rm -rf lib/cmake/cblas* lib/cmake/lapack* lib/pkgconfig
rm -rf lib/cmake lib/pkgconfig