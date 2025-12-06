#!/usr/bin/env bash
set -e

# detect OS
OS_NAME="$(uname -s)"
case "$OS_NAME" in
    Linux*)  OS="linux-gnu";;
    Darwin*) OS="apple-darwin";;
    Windows*) OS="w64-mingw32";;
    *) echo "Unsupported OS: $OS_NAME"; exit 1 ;;
esac
# detect architecture
ARCH_NAME="$(uname -m)"
case "$ARCH_NAME" in
    x86_64|amd64) ARCH="x86_64";;
    arm64|aarch64) ARCH="aarch64";;
    *) echo "Unknown architecture '$ARCH_NAME'."; exit 1;;
esac

# create and change directory
mkdir -p interfaces/Python/third_party
cd interfaces/Python/third_party

# download BQPD
VERSION="v1.0.0"
ASSET_URL="https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-${VERSION}%2B0/BQPD.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
echo "Downloading: ${ASSET_URL}"
mkdir bqpd && cd bqpd
curl -L -o bqpd.tar.gz "$ASSET_URL"
tar xf bqpd.tar.gz
ls -l
pwd
cd ..

# download MUMPS
VERSION="v5.8.0"
ASSET_URL="https://github.com/amontoison/MUMPS_static_jll.jl/releases/download/MUMPS_static-${VERSION}%2B0/MUMPS_static.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
echo "Downloading: ${ASSET_URL}"
mkdir mumps && cd mumps
curl -L -o mumps.tar.gz "$ASSET_URL"
tar xf mumps.tar.gz
ls -l
pwd
cd ..