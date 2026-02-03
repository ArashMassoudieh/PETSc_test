cd ~/Projects/PETSc_test

cat > diag_std_headers.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

echo "== PWD =="
pwd
echo

echo "== mpicxx =="
command -v mpicxx || true
mpicxx -show || true
echo

echo "== possible shadowing files (string/cstddef) in repo =="
find . -maxdepth 4 \( -name string -o -name string.h -o -name cstddef -o -name cstddef.h \) -print || true
echo

echo "== include trace for <string> (first 60 lines) =="
echo '#include <string>' | mpicxx -x c++ -std=c++17 -E -H - 2>&1 | head -n 60
echo

echo "== compile sanity test (std::string) =="
cat <<'CPP' | mpicxx -x c++ -std=c++17 -o /tmp/stdtest /dev/stdin && echo "Sanity compile OK"
/* sanity */
#include <string>
int main(){ std::string s="ok"; return (int)s.size(); }
CPP
EOF

chmod +x diag_std_headers.sh
./diag_std_headers.sh

