#!/usr/bin/env bash
set -euo pipefail

PETSC_DIR="${PETSC_DIR:-/home/arash/Projects/petsc}"
PETSC_ARCH="${PETSC_ARCH:-arch-linux-c-opt}"

LIBDIR="$PETSC_DIR/$PETSC_ARCH/lib"

echo "=== PETSc dynamic-link check (ldd) ==="
echo "PETSC_DIR = $PETSC_DIR"
echo "PETSC_ARCH= $PETSC_ARCH"
echo "LIBDIR    = $LIBDIR"
echo

# Find a PETSc shared library to inspect
CANDIDATES=()
[[ -f "$LIBDIR/libpetsc.so" ]] && CANDIDATES+=("$LIBDIR/libpetsc.so")
CANDIDATES+=( "$LIBDIR"/libpetsc.so.* )

TARGET=""
for f in "${CANDIDATES[@]}"; do
  if [[ -f "$f" ]]; then TARGET="$f"; break; fi
done

if [[ -z "$TARGET" ]]; then
  echo "ERROR: Could not find libpetsc.so in:"
  echo "  $LIBDIR"
  echo "Maybe PETSc was built static-only."
  echo "If static, use the petscvariables check (script #1) and/or rebuild with shared libs."
  exit 2
fi

echo "Inspecting: $TARGET"
echo

echo "--- ldd (filtered for blas/openblas/mkl/blis/lapack) ---"
ldd "$TARGET" | grep -iE 'openblas|blas|lapack|mkl|blis|atlas|accelerate' || true
echo

echo "--- full ldd (last 80 lines) ---"
ldd "$TARGET" | tail -n 80 || true
echo

if ldd "$TARGET" | grep -qi openblas; then
  echo "✅ Definitive: libpetsc.so links to OpenBLAS."
else
  echo "❌ Definitive: no OpenBLAS in libpetsc.so dependencies."
  echo "   (It might be using system BLAS/LAPACK, MKL, or static BLAS.)"
fi

