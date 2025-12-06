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
pwd

# download BQPD
VERSION="v1.0.0"
REPO="https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-${VERSION}%2B0"
ASSET_NAME="BQPD.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
mkdir bqpd
curl -L -o bqpd/bqpd.tar.gz "$ASSET_URL"
tar -xzvf bqpd/bqpd.tar.gz

# download MUMPS
VERSION="v5.8.0"
REPO="https://github.com/amontoison/MUMPS_static_jll.jl/releases/download/MUMPS_static-${VERSION}%2B0"
ASSET_NAME="MUMPS_static.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
mkdir mumps
curl -L -o mumps/mumps.tar.gz "$ASSET_URL"
tar -xzvf mumps/mumps.tar.gz

# download HiGHS
VERSION="v1.11.0"
REPO="https://github.com/amontoison/HiGHS_static_jll.jl/releases/download/HiGHS_static-${VERSION}%2B0"
ASSET_NAME="HiGHS_static.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"
ASSET_URL="${REPO}/${ASSET_NAME}"
echo "Downloading: ${ASSET_URL}"
mkdir highs
curl -L -o highs/highs.tar.gz "$ASSET_URL"
tar -xzvf highs/highs.tar.gz