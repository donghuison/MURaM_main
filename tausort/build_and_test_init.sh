#!/usr/bin/env bash
set -euo pipefail

DIR="/Users/donghuison/workspace/myGit/MURaM-study/MURaM_main/tausort"
cd "$DIR"

# Ensure required data files exist (Nlam=328 path)
for f in G2_1D.dat p00big2_asplund.bdf asplund_abs_cont.dat asplund_sca_cont.dat; do
  if [[ ! -f "$f" ]]; then
    echo "Missing required data: $f" >&2
    exit 1
  fi
done

# Build C objects without conflicting main
gcc -O3 -std=c11 -Dmain=tausort_main -c "$DIR/tausort.c" -o "$DIR/tausort.o"

# Build Fortran driver and link
gfortran -O3 -std=f2008 "$DIR/test_init_driver.f90" "$DIR/tausort.o" -o "$DIR/run_tausort_init" -lm

# Run and capture output
./run_tausort_init | tee out.txt

# Verify expected banners around initialize()
grep -q "1    Input Finished. Beginning initialisation." out.txt
grep -q "Initialisation complete." out.txt
echo "[OK] Fortran driver successfully called C input()/initialize()."


