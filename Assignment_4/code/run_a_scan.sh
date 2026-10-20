#!/bin/bash
# Part A, time-step explanation: scan of the impact speed at fixed time step.
# The exact e of the linear model does not depend on the impact speed, so any
# change of e with v at fixed dt is a property of the integration (see the
# notebook). 2 tests x 4 DT_FRACTION x 31 speeds (0.10 to 1.00 m/s, step 0.03)
# = 248 runs, in parallel, one directory each: runs/A_scan_<test>_f<frac>_v<v>/.
# Run from Assignment_4/:   bash code/run_a_scan.sh      (analysed by code/analyse_a.py)
# On Windows/MSYS2 first:   export PATH=/c/msys64/ucrt64/bin:$PATH
one_run() {
  test=$1; frac=$2; v=$3
  d=runs/A_scan_${test}_f${frac}_v${v}
  mkdir -p $d
  gcc -O3 -DRUN=1 -DNUM_TEST=$test -DV_IMPACT=$v -DDT_FRACTION=$frac code/*.c -o $d/dem.exe -lm
  (cd $d && ./dem.exe > run.dat)
}
export -f one_run
SECONDS=0
for test in 1 2; do
  for frac in 0.2 0.1 0.05 0.025; do
    for i in $(seq 10 3 100); do
      printf "%s %s %s\n" $test $frac $(awk -v i=$i 'BEGIN{printf "%.2f", i/100}')
    done
  done
done | xargs -P 16 -L 1 bash -c 'one_run $0 $1 $2'
echo "248 scan runs, total wall-clock (compile + run, 16 in parallel): ${SECONDS} s"
