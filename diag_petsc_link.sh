#!/usr/bin/env bash
set -euo pipefail

echo "=== PETSc link diagnostic ==="
echo "PWD: $(pwd)"
echo

# 1) Basic build context
echo "== Git / repo root =="
git rev-parse --show-toplevel 2>/dev/null || true
echo

echo "== Qt / qmake =="
command -v qmake >/dev/null && qmake -v || echo "qmake not found"
echo

echo "== MPI compiler wrapper =="
command -v mpicxx >/dev/null && which mpicxx || echo "mpicxx not found"
echo

# 2) PETSc env + expected paths
echo "== PETSc env =="
echo "PETSC_DIR=$PETSC_DIR"
echo "PETSC_ARCH=$PETSC_ARCH"
echo

# 3) Find petscsys.h (proves PETSc include root)
echo "== Find petscsys.h =="
find /home/arash/Projects -maxdepth 4 -name petscsys.h 2>/dev/null | head -n 20 || true
echo

# 4) Find libpetsc and check it contains PetscInitialized
echo "== Find PETSc libraries (libpetsc*) =="
PETSCLIBS=$(find /home/arash/Projects -maxdepth 6 -type f \( -name "libpetsc.so" -o -name "libpetsc.so.*" -o -name "libpetsc.a" \) 2>/dev/null || true)
if [[ -z "$PETSCLIBS" ]]; then
  echo "No libpetsc found under /home/arash/Projects. PETSc might not be built/installed where you think."
else
  echo "$PETSCLIBS" | head -n 20
fi
echo

# Choose one libpetsc for symbol tests (prefer .so)
PETSCLIB_CHOSEN=""
if [[ -n "$PETSCLIBS" ]]; then
  PETSCLIB_CHOSEN=$(echo "$PETSCLIBS" | grep -E '\.so(\.|$)' | head -n 1 || true)
  [[ -z "$PETSCLIB_CHOSEN" ]] && PETSCLIB_CHOSEN=$(echo "$PETSCLIBS" | head -n 1)
fi

echo "== Check PetscInitialized symbol inside libpetsc =="
if [[ -z "$PETSCLIB_CHOSEN" ]]; then
  echo "SKIP: no libpetsc chosen"
else
  echo "Using: $PETSCLIB_CHOSEN"
  if file "$PETSCLIB_CHOSEN" | grep -qi "shared object"; then
    # dynamic symbol table
    nm -D "$PETSCLIB_CHOSEN" 2>/dev/null | grep -E "PetscInitialized$" && echo "FOUND" || echo "NOT FOUND"
  else
    # static archive
    nm "$PETSCLIB_CHOSEN" 2>/dev/null | grep -E "PetscInitialized$" && echo "FOUND" || echo "NOT FOUND"
  fi
fi
echo

# 5) If petscvariables exists, print key vars (this is what qmake includes)
echo "== petscvariables check =="
# Try to infer PETSC_DIR/ARCH if env not set (Jason default)
PETSC_DIR_GUESS="${PETSC_DIR:-/home/arash/Projects/petsc}"
PETSC_ARCH_GUESS="${PETSC_ARCH:-arch-linux-c-opt}"
PV="$PETSC_DIR_GUESS/$PETSC_ARCH_GUESS/lib/petsc/conf/petscvariables"
echo "Guess petscvariables: $PV"
if [[ -f "$PV" ]]; then
  echo "Found petscvariables."
  echo "-- grep PETSC_LIB / PETSC_EXTERNAL_LIB* --"
  grep -E '^(PETSC_LIB|PETSC_EXTERNAL_LIB|PETSC_EXTERNAL_LIB_BASIC) *=' "$PV" || true
else
  echo "NOT FOUND: $PV"
fi
echo

# 6) Look for Qt shadow-build confusion (your error path shows build-PETCs_test_Behzad...)
echo "== Look for build directories with old names (Behzad etc.) =="
ls -d ../build-* 2>/dev/null | head -n 50 || true
echo

# 7) If a Makefile exists in current dir, print the final link command used
echo "== Try to extract link command (make -n / make V=1) =="
if [[ -f "Makefile" ]]; then
  echo "-- make -n (show commands, do not run) --"
  make -n 2>/dev/null | grep -E "(mpicxx|g\+\+|c\+\+).*-o .*" | tail -n 5 || true
  echo
  echo "-- make V=1 (runs; may fail; shows full link line) --"
  set +e
  make clean >/dev/null 2>&1
  make -j1 V=1 2>&1 | tee /tmp/build_v1.log
  set -e
  echo
  echo "-- Extract likely link lines from /tmp/build_v1.log --"
  grep -E "(mpicxx|g\+\+|c\+\+).*-o .*" /tmp/build_v1.log | tail -n 10 || true
else
  echo "No Makefile here. Run this script inside your Qt build directory OR run qmake first."
fi
echo

echo "=== Interpretation tips ==="
cat <<'EOF'
1) If "Check PetscInitialized symbol" says NOT FOUND:
   -> You are linking the wrong PETSc library (or PETSc build is broken).

2) If it says FOUND but the link still fails:
   -> PETSc libs are NOT on the final link line OR appear before objects (order issue).
   -> Check the extracted link line: it must contain -lpetsc (or the PETSC_LIB expansion).

3) If the build directory list shows something like build-PETCs_test_Behzad-Desktop-Debug:
   -> Qt Creator is compiling in an old shadow-build folder. Create a new build folder for Jason and rerun qmake.

4) If petscvariables exists but grep shows empty PETSC_LIB:
   -> petscvariables not generated correctly or not the correct one for your build.
EOF

