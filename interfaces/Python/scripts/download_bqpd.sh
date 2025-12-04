#!/usr/bin/env bash
set -e

mkdir -p third_party/bqpd
cd third_party/bqpd

VERSION="v1.0.0"

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

ASSET_URL="https://github.com/leyffer/BQPD_jll.jl/releases/download/BQPD-${VERSION}%2B0/BQPD.${VERSION}.${ARCH}-${OS}-libgfortran5.tar.gz"

echo "Downloading: ${ASSET_URL}"
curl -L -o bqpd.tar.gz "$ASSET_URL"
tar xf bqpd.tar.gz

ls -l
