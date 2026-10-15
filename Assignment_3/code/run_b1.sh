#!/bin/bash
# B1: energy conservation with the conservative force only (dissipative and random force off).
# Run from the Assignment_3 directory with gcc on the PATH:   bash code/run_b1.sh
#
# Start: the end state of the thermostatted monomer fluid of C1 (restart_files/C1_a25.dat, all forces on,
# k_BT = 1), so the run starts from the state that is simulated everywhere else. Every run lasts t = 8.
# Output: data/B1_<tag>.csv with columns step,time,Epot,Ekin,Etot.
set -e

START=restart_files/C1_a25.dat
[ -f $START ] || { echo "missing $START: run  RUNS=C1_a25 bash code/run_c.sh  first"; exit 1; }

for entry in 004:0.04 002:0.02 001:0.01 0005:0.005; do
    tag=${entry%%:*}; dt=${entry##*:}
    dir=build_b1/$tag; rm -rf $dir; mkdir -p $dir
    cp code/*.c code/*.h $dir/
    steps=$(python -c "print(round(8/$dt))")
    python - $dir/setparameters.c $dt $steps $START <<'PY'
import re, sys
f, dt, steps, start = sys.argv[1:]
s = open(f, encoding='utf8', newline='').read()
def setp(k, v):
    global s
    s, n = re.subn(r'(p_parameters->%s\s*=\s*)[^;]*;' % re.escape(k), lambda m: m.group(1) + v + ';', s, count=1)
    assert n == 1, k
for k, v in [('dissipative_random_force_on', '0'), ('conservative_force_on', '1'), ('chain_length', '1'),
             ('a[1]', '25'), ('demix_initialisation', '0'), ('dt', dt), ('num_dt_steps', steps),
             ('load_restart', '1'), ('r_shell', '0.4'), ('num_dt_pdb', '100000000'), ('num_dt_restart', '100000000'),
             ('composition_profile_on', '0'), ('profile_on', '0'), ('analysis_on', '0')]:
    setp(k, v)
s, n = re.subn(r'p_parameters->L\s*=\s*\(struct Vec3D\)\{[^}]*\};', 'p_parameters->L = (struct Vec3D){10, 10, 10};', s, count=1)
assert n == 1
for key, val in (('restart_in_filename', start), ('restart_out_filename', 'restart_files/B1_unused.dat'), ('filename_pdb', 'build_b1/unused')):
    pat = re.compile(r'(?:strcpy|snprintf)\(\s*p_parameters->%s\s*,.*?\);' % key, re.S)
    s, n = pat.subn(lambda m: 'strcpy(p_parameters->%s, "%s");' % (key, val), s, count=1)
    assert n == 1, key
open(f, 'w', encoding='utf8', newline='').write(s)
PY
    gcc -O3 $dir/*.c -o $dir/dpd -lm
    ./$dir/dpd > data/B1_$tag.csv
    echo "B1_$tag.csv written (dt = $dt, $steps steps)"
done
rm -rf build_b1 restart_files/B1_unused.dat build_b1
