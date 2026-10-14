#!/bin/bash
# Part C production runs (C1 g(r), C2/C3/C4 demixing and surface tension).
# Run from the Assignment_3 directory with gcc on the PATH:   bash code/run_c.sh
#
# All run settings live in setparameters.c, so every run is built from a patched
# copy of the sources in build_c/<run>/ (code/setparameters.c itself is untouched).
# Output goes to data/<run>_*.csv, restart files to restart_files/<run>.dat.
set -e

RUNS=${RUNS:-"C1_a25 C2_a37 C3_a35 C3_a39 C3_a43 C3_a47 C4_a25"}
STEPS=${STEPS:-60000}
BLOCK=${BLOCK:-5000}   # steps per profile/g(r) block

patch_param () {   # patch_param <file> key=value ...   (s:key=text for strcpy strings)
    python - "$@" <<'PY'
import re, sys
f = sys.argv[1]
s = open(f, encoding='utf8', newline='').read()
for a in sys.argv[2:]:
    k, v = a.split('=', 1)
    if k.startswith('s:'):
        k = k[2:]
        pat = re.compile(r'(strcpy\(\s*p_parameters->%s,\s*)"[^"]*"' % re.escape(k))
        s, n = pat.subn(lambda m: m.group(1) + '"' + v + '"', s, count=1)
    elif k == 'L':   # box: "Lx,Ly,Lz"
        pat = re.compile(r'p_parameters->L\s*=\s*\(struct Vec3D\)\{[^}]*\};')
        s, n = pat.subn(lambda m: 'p_parameters->L = (struct Vec3D){%s};' % v, s, count=1)
    elif k == 'filename_analysis':
        pat = re.compile(r'snprintf\(\s*p_parameters->filename_analysis.*?\);', re.S)
        s, n = pat.subn(lambda m: 'strcpy(p_parameters->filename_analysis, "%s");' % v, s, count=1)
    else:
        pat = re.compile(r'(p_parameters->%s\s*=\s*)[^;]*;' % re.escape(k))
        s, n = pat.subn(lambda m: m.group(1) + v + ';', s, count=1)
    assert n == 1, k
open(f, 'w', encoding='utf8', newline='').write(s)
PY
}

build () {   # build <run> <a_AB> <demix> <box> <steps> <extra patches...>
    local run=$1 aab=$2 demix=$3 box=$4 steps=$5; shift 5
    local dir=build_c/$run
    mkdir -p $dir
    cp code/*.c code/*.h $dir/
    patch_param $dir/setparameters.c "a[1]=$aab" demix_initialisation=$demix "L=$box" \
        load_restart=0 r_shell=0.4 chain_length=1 conservative_force_on=1 dissipative_random_force_on=1 \
        composition_profile_on=0 num_dt_steps=$steps num_dt_pdb=$steps num_dt_restart=$steps \
        s:filename_pdb=data/${run}_traj s:restart_out_filename=restart_files/${run}.dat \
        filename_analysis=data/${run} "$@"
    gcc -O3 $dir/*.c -o $dir/dpd -lm
}

mkdir -p build_c data restart_files

for run in $RUNS; do
    aab=${run##*_a}
    case $run in
        C1_*)   # monomer fluid: cubic box, g(r) after 4000 steps of equilibration, blocks of 4000 steps
            build $run $aab 0 "10,10,10" 20000 analysis_on=1 num_dt_phi=20 analysis_start_step=4000 num_dt_block=4000 ;;
        C2_*)   # Fig. 6: snapshots every 10000 steps to follow the demixing
            build $run $aab 1 "20,10,10" $STEPS profile_on=1 num_dt_block=$BLOCK num_dt_pdb=10000 ;;
        *)      # C3/C4 series and the control
            build $run $aab 1 "20,10,10" $STEPS profile_on=1 num_dt_block=$BLOCK ;;
    esac
    # (num_part = rho * Lx*Ly*Lz is computed in setparameters.c from L)
done

for run in $RUNS; do
    ( start=$(date +%s)
      ./build_c/$run/dpd | awk 'NR==1 || NR%50==0' > data/${run}_energy.csv
      echo "$run finished in $(( $(date +%s) - start )) s" ) &
done
wait
echo "all runs finished"
