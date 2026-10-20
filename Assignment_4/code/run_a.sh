#!/bin/bash
# Part A: the 24 restitution runs (NUM_TEST x V_IMPACT x DT_FRACTION).
# Run from Assignment_4/:   bash code/run_a.sh
# One directory per run, runs/A_<test>_v<v>_f<frac>/, holding dem.exe and run.dat
# (the table printed every step). On Windows/MSYS2 first:
#   export PATH=/c/msys64/ucrt64/bin:$PATH
set -e
SECONDS=0
for test in 1 2; do
  for v in 0.1 0.3 1.0; do
    for frac in 0.2 0.1 0.05 0.025; do
      d=runs/A_${test}_v${v}_f${frac}
      mkdir -p $d
      gcc -O3 -Wall -Wextra -DRUN=1 -DNUM_TEST=$test -DV_IMPACT=$v -DDT_FRACTION=$frac code/*.c -o $d/dem.exe -lm
      (cd $d && ./dem.exe > run.dat)
    done
  done
done
echo "24 runs, total wall-clock (compile + run): ${SECONDS} s"
