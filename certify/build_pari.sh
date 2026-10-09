#!/usr/bin/env bash
# Build PARI/GP 2.19.0 with the bnftestprimes_range patch into ./.tools/pari.
# The patch adds a range-restricted copy of PARI's own bnftestprimes() (Phase 1
# of bnfcertify) so the prime loop can be split across processes.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"
VER=2.19.0
SHA256=$(cat certify/pari-$VER.sha256)
PREFIX="$ROOT/.tools/pari"
SRC="$ROOT/.tools/src"

if [[ -x "$PREFIX/bin/gp" ]]; then echo "already built: $PREFIX/bin/gp"; exit 0; fi
mkdir -p "$SRC"
cd "$SRC"
if [[ ! -f pari-$VER.tar.gz ]]; then
  curl -fsSL -o pari-$VER.tar.gz "https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-$VER.tar.gz"
fi
echo "$SHA256  pari-$VER.tar.gz" | shasum -a 256 -c -
rm -rf pari-$VER && tar xzf pari-$VER.tar.gz
cd pari-$VER
patch -p1 < "$ROOT/certify/pari-$VER-bnftestprimes-range.patch"

GMP_FLAG=()
if command -v brew >/dev/null && brew --prefix gmp >/dev/null 2>&1; then
  GMP_FLAG=(--with-gmp="$(brew --prefix gmp)")
fi
./Configure --prefix="$PREFIX" "${GMP_FLAG[@]}" > "$SRC/configure.out" 2>&1
make -j"$(sysctl -n hw.ncpu 2>/dev/null || nproc)" gp > "$SRC/make.out" 2>&1
make install > "$SRC/install.out" 2>&1
echo "built: $PREFIX/bin/gp"
