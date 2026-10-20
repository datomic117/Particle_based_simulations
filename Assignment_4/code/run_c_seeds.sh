#!/bin/bash
# Part C follow-up, job 1: the silo run (D = 7 d_p, defaults) for several random seeds, each with the edge of the hole ON
# (the shipped code) and OFF (the early `return false` of run_c_noedge.sh in a scratch copy). code/ itself is not changed.
# Seed 12345 is the seed of the shipped code (data/D7.dat, data/C_noedge_D7.dat) and is not repeated here.
# Run from Assignment_4/ with gcc on the PATH:   bash code/run_c_seeds.sh        (8 runs in parallel, about 5 min)
# Each run gets runs/C_seed<s>_<edge|noedge>/ (gitignored) with a scratch copy of the sources in src/ (constants.h with the new
# SEED, and for `noedge` walls.c with the edge off); the wall-clock goes to time.txt and the table is copied afterwards:
#   cp runs/C_seed<s>_edge/run.dat data/C_seed<s>_edge.dat     (same for noedge)
SEEDS="${SEEDS:-1 2 3 4}"
run() {  # seed  edge|noedge
  d=runs/C_seed$1_$2
  mkdir -p $d/src
  cp code/*.c code/*.h $d/src/
  sed -i "s|^#define SEED .*|#define SEED $1u|" $d/src/constants.h
  grep -c "^#define SEED $1u" $d/src/constants.h > /dev/null || { echo "seed not set for $d"; exit 1; }
  if [ "$2" = "noedge" ]; then
    sed -i 's|^    double ex = 1.0, ey = 0.0; // e_r; the arbitrary direction on the axis|    return false; /* edge switched off for the comparison */\n&|' $d/src/walls.c
    grep -c "edge switched off" $d/src/walls.c > /dev/null || { echo "edge not switched off in $d"; exit 1; }
  fi
  gcc -O3 -DRUN=3 $d/src/*.c -o $d/dem.exe -lm
  ( cd $d && { time ./dem.exe > run.dat; } 2> time.txt )
}
for s in $SEEDS; do
  run $s edge &
  run $s noedge &
done
wait
