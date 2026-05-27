#!/usr/bin/env bash
set -e

# detect OS
OS_NAME="$(uname -s)"
case "$OS_NAME" in
    Linux*)
      if ldd --version 2>&1 | grep -qi musl; then
        echo "Unsupported OS: linux-musl"
        exit 1
      else
        OS="linux-gnu"
        OS_KRYLOV="linux"
        EXTENSION_KRYLOV="tar.gz"
      fi
      ;;
    Darwin*)
		OS="apple-darwin"
		OS_KRYLOV="macos"
		EXTENSION_KRYLOV="tar.gz"
		;;
    Windows*|MINGW64_NT*|MSYS_NT*)
		OS="w64-mingw32"
		OS_KRYLOV="windows"
		EXTENSION_KRYLOV="zip"
		;;
    *)
		echo "Unsupported OS: $OS_NAME"
		exit 1
		;;
esac
# detect architecture
ARCH_NAME="$(uname -m)"
case "$ARCH_NAME" in
    x86_64|amd64)
		ARCH="x86_64"
		ARCH_KRYLOV="x86_64"
		;;
    arm64|aarch64)
		ARCH="aarch64";
		case "$OS_KRYLOV" in
			linux*)
				ARCH_KRYLOV="aarch64"
				;;
			macos*)
				ARCH_KRYLOV="arm64"
				;;
			*)
				echo "Unsupported aarch64 OS: $OS_KRYLOV"
				exit 1
				;;
		esac
		;;
    *)
		echo "Unknown architecture '$ARCH_NAME'."
		exit 1
		;;
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

# on MinGW, compile HiGHS on the fly and replace that of UnoUtils
if [[ "$OS" == "w64-mingw32" && "${UNO_TOOLCHAIN:-mingw}" == "mingw" ]]; then
	# delete UnoUtils' HiGHS
	rm -rf lib/libhighs* include/highs include/Highs*.h bin/highs*
	
	# The toolchain is C:/mingw64 via the pinned CMAKE_*_COMPILER/CC/CXX;
	# /mingw64/bin supplies only cmake/make
	export PATH="$PATH:/mingw64/bin"
	if ! command -v cmake >/dev/null 2>&1; then
		pacman -Sy --needed --noconfirm mingw-w64-x86_64-cmake
	fi

	# fetch HiGHS source
	BUILD_ROOT="$(mktemp -d)"
	mkdir "${BUILD_ROOT}/HiGHS"
	VERSION="v1.15.0"
	ASSET_URL="https://github.com/ERGO-Code/HiGHS/archive/refs/tags/${VERSION}.tar.gz"
	echo "Downloading: ${ASSET_URL}"
	curl -fL -o "${BUILD_ROOT}/HiGHS-src.tar.gz" "$ASSET_URL"
	tar -xzf "${BUILD_ROOT}/HiGHS-src.tar.gz" -C "${BUILD_ROOT}/HiGHS" --strip-components=1

	# Build HiGHS with the SAME compiler the consuming workflow links with.
	# A mismatch here is what produces the __emutls_v._ZSt*__once_call* and
	# libstdc++ undefined references (e.g. GCC 8.1.0 emutls vs GCC 16 native TLS).
	towin() { if command -v cygpath >/dev/null 2>&1; then cygpath -m "$1"; else printf '%s' "$1"; fi; }
	CC_BIN="${CMAKE_C_COMPILER:-${CC:-$(command -v gcc || true)}}"
	CXX_BIN="${CMAKE_CXX_COMPILER:-${CXX:-$(command -v g++ || true)}}"
	MAKE_BIN="${CMAKE_MAKE_PROGRAM:-$(command -v mingw32-make || command -v make || true)}"
	if [[ -z "$CC_BIN" || -z "$CXX_BIN" || -z "$MAKE_BIN" ]]; then
		echo "No MinGW toolchain found; set CMAKE_{C,CXX}_COMPILER/CMAKE_MAKE_PROGRAM or put gcc/g++/make on PATH"; exit 1
	fi
	echo "Building HiGHS with $CXX_BIN ($("$CXX_BIN" -dumpversion))"

	GEN_FLAGS=(-G "Unix Makefiles"
		-DCMAKE_MAKE_PROGRAM="$(towin "$MAKE_BIN")"
		-DCMAKE_C_COMPILER="$(towin "$CC_BIN")"
		-DCMAKE_CXX_COMPILER="$(towin "$CXX_BIN")")
		
	# native cmake needs Windows paths, not MSYS ones
	SRC_W="$(cygpath -m "${BUILD_ROOT}/HiGHS")"
	BLD_W="$(cygpath -m "${BUILD_ROOT}/build")"
	PREFIX_W="$(cygpath -m "${BUILD_ROOT}/install")"
	BLAS_W="$(cygpath -m "${PWD}/lib/libblas.a")"

	# build
	cmake -S "$SRC_W" -B "$BLD_W" \
		"${GEN_FLAGS[@]}" \
		-DCMAKE_INSTALL_PREFIX="$PREFIX_W" \
		-DCMAKE_BUILD_TYPE=Release \
		-DBUILD_SHARED_LIBS=OFF \
		-DZLIB=OFF \
		-DHIPO=ON \
		-DBUILD_EXAMPLES=OFF \
		-DBUILD_TESTING=OFF \
		-DBUILD_CXX_EXE=OFF \
		-DBLAS_LIBRARIES="$BLAS_W" \
		-DBUILD_SHARED_EXTRAS_LIB=OFF \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON

	cmake --build "$BLD_W" --config Release --parallel
	cmake --install "$BLD_W" --config Release

	cp -a "${BUILD_ROOT}/install/lib/."     lib
	cp -a "${BUILD_ROOT}/install/include/." include
	rm -rf "${BUILD_ROOT}"
fi

# download Krylov.jl
VERSION="v0.0.23"
REPO="https://github.com/amontoison/Krylov.jl/releases/download/${VERSION}"
ASSET_NAME="libkrylov-${OS_KRYLOV}-${ARCH_KRYLOV}.${EXTENSION_KRYLOV}"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
if [[ "$OS_KRYLOV" == "windows" ]]; then
	curl -L -o libkrylov.zip "$ASSET_URL"
	unzip libkrylov.zip
else
	curl -L -o libkrylov.tar.gz "$ASSET_URL"
	tar -xzf libkrylov.tar.gz
fi
pwd

# delete unwanted directories
# rm -rf lib/cmake/cblas* lib/cmake/lapack* lib/pkgconfig
rm -rf lib/cmake lib/pkgconfig