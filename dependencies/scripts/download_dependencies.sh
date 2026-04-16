#!/usr/bin/env bash
set -e

# detect OS
OS_NAME="$(uname -s)"
case "$OS_NAME" in
    Linux*)  OS="linux-gnu";;
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
VERSION="1.0.0"
REPO="https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-v${VERSION}%2B0"
ASSET_NAME="BQPD.v${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
curl -L -o BQPD.tar.gz "$ASSET_URL"
tar -xzvf BQPD.tar.gz
pwd

# download UnoUtils: MUMPS (+ METIS, BLAS and LAPACK) and HiGHS
VERSION="2026.4.16"
REPO="https://github.com/amontoison/UnoUtils_jll.jl/releases/download/UnoUtils-v${VERSION}%2B0"
ASSET_NAME="UnoUtils.v${VERSION}.${ARCH}-${OS}-libgfortran5-cxx11.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
curl -L -o UnoUtils.tar.gz "$ASSET_URL"
tar -xzvf UnoUtils.tar.gz
pwd

# delete unwanted directories
# rm -rf lib/cmake/cblas* lib/cmake/lapack* lib/pkgconfig
rm -rf lib/cmake lib/pkgconfig
