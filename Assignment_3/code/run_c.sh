#!/bin/bash
# C1/C2 production runs: NREP independent replicas, each an NVT equilibration
# (Berendsen, tau = 1 ps) followed by an NVE production run, all in parallel.
# Run from this directory (code/) with gcc on the PATH:   bash run_c.sh
#
# Run configuration lives in setparameters.c only, so every binary is built
# from a patched copy of it (setparameters.c itself is restored afterwards).
set -e

NREP=${NREP:-5}
EQUIL_STEPS=${EQUIL_STEPS:-200000}   # 200 ps at 1 fs
PROD_STEPS=${PROD_STEPS:-400000}     # 400 ps at 1 fs per replica
TAU_1PS=0.9118                       # 1 ps in internal time units

cp setparameters.c setparameters.c.orig

patch_param () {   # patch_param key value   |   patch_param s:key text
    python - "$@" <<'PY'
import re, sys
s = open('setparameters.c', encoding='utf8', newline='').read()
for a in sys.argv[1:]:
    k, v = a.split('=', 1)
    if k.startswith('s:'):
        k = k[2:]
        pat = re.compile(r'(strcpy\(\s*p_parameters->%s,\s*)"[^"]*"' % k)
        s, n = pat.subn(lambda m: m.group(1) + '"' + v + '"', s, count=1)
    else:
        pat = re.compile(r'(p_parameters->%s\s*=\s*)[^;]*;' % k)
        s, n = pat.subn(lambda m: m.group(1) + v + ';', s, count=1)
    assert n == 1, k
open('setparameters.c', 'w', encoding='utf8', newline='').write(s)
PY
}

for k in $(seq 1 $NREP); do
    # Equilibration: new velocities per replica, thermostat on, dihedral statistics in 10 ps blocks
    cp setparameters.c.orig setparameters.c
    patch_param load_restart=1 reseed_velocities=$((100 + k)) is_NVT=1 tau_T=$TAU_1PS \
        num_dt_steps=$EQUIL_STEPS num_dt_output=100 num_dt_pdb=100000000 num_dt_restart=$EQUIL_STEPS \
        analysis_on=1 msd_on=0 num_dt_block=10000 \
        s:filename_pdb=c_pdb_equil_$k s:restart_in_filename=restart_files/restart_b9_tau001.dat \
        s:restart_out_filename=restart_files/restart_c_equil_$k.dat s:filename_analysis=../data/c_equil_$k
    gcc -O3 *.c -o md_equil_$k -lm

    # Production: NVE from the equilibrated state, dihedral statistics and MSD in 100 ps blocks
    cp setparameters.c.orig setparameters.c
    patch_param load_restart=1 reseed_velocities=0 is_NVT=0 \
        num_dt_steps=$PROD_STEPS num_dt_output=100 num_dt_pdb=100000000 num_dt_restart=$PROD_STEPS \
        analysis_on=1 msd_on=1 num_dt_block=100000 \
        s:filename_pdb=c_pdb_prod_$k s:restart_in_filename=restart_files/restart_c_equil_$k.dat \
        s:restart_out_filename=restart_files/restart_c_prod_$k.dat s:filename_analysis=../data/c_prod_$k
    gcc -O3 *.c -o md_prod_$k -lm
done
mv setparameters.c.orig setparameters.c

for k in $(seq 1 $NREP); do
    ( ./md_equil_$k > ../data/c_equil_$k.csv && ./md_prod_$k > ../data/c_prod_$k.csv ) &
done
wait
rm -f md_equil_* md_prod_* c_pdb_*.pdb
echo "all replicas finished"
