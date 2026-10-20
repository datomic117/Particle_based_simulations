#!/bin/bash
# Part C2 comparison: the silo run (D = 7 d_p, defaults) with the edge of the hole switched off, i.e. as the starting code behaves
# before Task C2 is done (a particle over the open hole never touches the plate). The edge branch of silo_bottom_wall in a scratch
# copy of the sources gets an early `return false`; code/ itself is not changed. Run from Assignment_4/ (gcc on the PATH):
#   bash code/run_c_noedge.sh
# Output: runs/C_noedge/run.dat (the table), copied to data/C_noedge_D7.dat by hand.
set -e
mkdir -p runs/C_noedge/src
cp code/*.c code/*.h runs/C_noedge/src/
sed -i 's|^    double ex = 1.0, ey = 0.0; // e_r; the arbitrary direction on the axis|    return false; /* edge switched off for the comparison */\n&|' runs/C_noedge/src/walls.c
grep -c "edge switched off" runs/C_noedge/src/walls.c
gcc -O3 -DRUN=3 runs/C_noedge/src/*.c -o runs/C_noedge/dem.exe -lm
(cd runs/C_noedge && time ./dem.exe > run.dat)
