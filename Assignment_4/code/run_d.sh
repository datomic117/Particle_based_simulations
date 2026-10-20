#!/bin/bash
# Part D3/D4: the extra silo runs (orifice D_ORIFICE in particle diameters, discharge time T_DISCHARGE in s).
# Run from Assignment_4/ with gcc on the PATH:  bash code/run_d.sh            (all four, in parallel)
# Each run gets its own directory runs/<name>/ (gitignored); the wall-clock goes to runs/<name>/time.txt and the table is
# copied to data/<name>.dat afterwards (cp runs/<name>/run.dat data/<name>.dat). D7 (the default) is made in Part C.
run() {  # name  D_ORIFICE  T_DISCHARGE
  mkdir -p runs/$1
  gcc -O3 -Wall -Wextra -DRUN=3 -DD_ORIFICE=$2 -DT_DISCHARGE=$3 code/*.c -o runs/$1/dem.exe -lm
  ( cd runs/$1 && { time ./dem.exe > run.dat; } 2> time.txt )
}
run D6 6 1.5 &
run D8 8 1.0 &
run D9 9 0.8 &
run D3 3 0.8 &
wait
